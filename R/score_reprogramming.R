#' Score drug candidates for beta-cell reprogramming potential
#'
#' Domain-specific scoring system for beta-cell reprogramming drug screening.
#' Designed as a generalizable framework: beta-cell reprogramming pathways
#' are used as defaults, but any cell type can be supported via custom
#' pathway dictionaries.
#'
#' Two scoring modes are provided for sensitivity analysis:
#' \itemize{
#'   \item \code{"conservative"}: connectivity (60\%) + validated pathways (40\%)
#'   \item \code{"comprehensive"}: connectivity (50\%) + validated pathways (30\%)
#'         + novel pathways (20\%)
#' }
#'
#' Connectivity score combines NCS magnitude (75\%) and WTCS_FDR significance
#' (25\%) for more robust ranking than NCS alone.
#'
#' @param forward_results Dataframe from query_cmap(). Must contain columns
#'   'NCS', 'WTCS_FDR', and an MOA column ('moa', 'MOA', 'moa_x', or 'moa.x').
#' @param mode Scoring mode: \code{"conservative"} (default) or
#'   \code{"comprehensive"}.
#' @param weight_ncs Numeric. Weight for connectivity score. Defaults: 0.6
#'   (conservative) or 0.5 (comprehensive).
#' @param weight_pathway Numeric. Weight for validated pathway score. Defaults:
#'   0.4 (conservative) or 0.3 (comprehensive).
#' @param weight_novel Numeric. Weight for novel pathway score (comprehensive
#'   mode only). Default: 0.2.
#' @param custom_validated_pathways Optional named list of pathway -> MOA
#'   keyword vectors. Overrides built-in beta-cell pathway dictionary.
#'   Use this to adapt ReprogDrugR to other cell types.
#' @param custom_novel_pathways Optional named list. Overrides built-in novel
#'   pathway dictionary.
#' @param toxicity_blacklist Character vector of additional drug names to flag
#'   as toxic (merged with built-in beta-cell toxicity list).
#' @param workers Integer. Number of parallel workers (default 1).
#'   Increase on servers (e.g., workers = 25).
#' @param reverse_only Logical. If TRUE (default), automatically retains only
#'   NCS < 0 compounds, which represent perturbations opposite to the query
#'   disease signature and are true disease-reversal candidates. NCS > 0
#'   compounds mimic the disease state and are excluded from scoring.
#'   If FALSE, all compounds are kept in the output, but the reprogramming
#'   score still only rewards NCS < 0 compounds — NCS > 0 compounds will
#'   receive score_ncs = 0 and score_fdr = 0 and thus score near zero.
#'
#' @return A dataframe sorted by reprogramming_score (descending) with columns:
#'   \itemize{
#'     \item \code{score_ncs}: Normalized absolute NCS for NCS < 0 compounds
#'       (0 to 1); NCS > 0 compounds receive 0
#'     \item \code{score_fdr}: FDR-based significance score (0 to 1);
#'       only awarded to NCS < 0 compounds
#'     \item \code{score_connectivity}: Combined NCS + FDR score (0 to 1)
#'     \item \code{score_pathway}: Validated pathway score (0 to 1)
#'     \item \code{matched_pathway}: Names of matched validated pathways
#'     \item \code{score_novel}: Novel pathway score (comprehensive only)
#'     \item \code{matched_novel}: Names of matched novel pathways
#'       (comprehensive only)
#'     \item \code{reprogramming_score}: Final weighted score
#'     \item \code{is_toxic}: Logical, TRUE if on beta-cell toxicity blacklist
#'     \item \code{scoring_mode}: Records which mode was used
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Recommended: pass full query_cmap() output, let reverse_only handle filtering
#' results_cons <- score_reprogramming(results_forward,
#'                                     mode = "conservative")
#'
#' # Comprehensive mode (sensitivity analysis)
#' results_comp <- score_reprogramming(results_forward,
#'                                     mode = "comprehensive")
#'
#' # Keep all compounds in output (NCS > 0 still scores near zero)
#' results_all  <- score_reprogramming(results_forward,
#'                                     mode = "conservative",
#'                                     reverse_only = FALSE)
#'
#' # Compare modes (Spearman correlation)
#' compare_modes(results_cons, results_comp)
#' }
score_reprogramming <- function(forward_results,
                                mode              = c("conservative",
                                                      "comprehensive"),
                                weight_ncs        = NULL,
                                weight_pathway    = NULL,
                                weight_novel      = NULL,
                                custom_validated_pathways = NULL,
                                custom_novel_pathways     = NULL,
                                toxicity_blacklist        = NULL,
                                workers                   = 1,
                                reverse_only              = TRUE) {

  mode <- match.arg(mode)

  # ---------- Input validation ----------
  if (!is.data.frame(forward_results) || nrow(forward_results) == 0)
    stop("forward_results must be a non-empty dataframe from query_cmap()")
  if (!("NCS" %in% colnames(forward_results)))
    stop("forward_results must contain an 'NCS' column")
  if (!is.numeric(workers) || workers < 1)
    stop("workers must be a positive integer")

  results <- forward_results

  # ---------- Direction filter ----------
  # In CMap/LINCS:
  #   NCS < 0 = perturbation opposite to query signature = disease-reversal candidate
  #   NCS > 0 = perturbation similar to query signature  = disease-mimicking compound
  # Only NCS < 0 compounds are valid reprogramming candidates.
  if (isTRUE(reverse_only)) {
    n_before  <- nrow(results)
    results   <- results[!is.na(results$NCS) & results$NCS < 0, , drop = FALSE]
    n_removed <- n_before - nrow(results)
    if (n_removed > 0) {
      warning(sprintf(
        paste0(
          "%d compound(s) with NCS > 0 removed. ",
          "These compounds mimic the disease signature and are not ",
          "reprogramming candidates. ",
          "Set reverse_only = FALSE to retain them in the output ",
          "(they will still score near zero)."
        ),
        n_removed
      ))
    }
    if (nrow(results) == 0)
      stop(paste0(
        "No NCS < 0 compounds found after direction filter. ",
        "Check that forward_results contains disease-reversal compounds. ",
        "Set reverse_only = FALSE to disable this filter."
      ))
  } else {
    # Keep all compounds but warn about NCS > 0.
    # Note: NCS > 0 compounds will receive score_ncs = 0 and score_fdr = 0
    # and will therefore score near zero. They are retained for reference only.
    n_positive <- sum(results$NCS > 0, na.rm = TRUE)
    if (n_positive > 0) {
      warning(sprintf(
        paste0(
          "%d compound(s) with NCS > 0 detected. ",
          "These compounds mimic the disease signature and will receive ",
          "score_ncs = 0 and score_fdr = 0. ",
          "They are retained in output for reference only."
        ),
        n_positive
      ))
    }
  }

  # ---------- Default weights ----------
  # Conservative:   connectivity 60% + validated pathway 40%
  # Comprehensive:  connectivity 50% + validated pathway 30% + novel 20%
  if (mode == "conservative") {
    if (is.null(weight_ncs))     weight_ncs     <- 0.60
    if (is.null(weight_pathway)) weight_pathway <- 0.40
    weight_novel <- 0
    total_w <- weight_ncs + weight_pathway
  } else {
    if (is.null(weight_ncs))     weight_ncs     <- 0.50
    if (is.null(weight_pathway)) weight_pathway <- 0.30
    if (is.null(weight_novel))   weight_novel   <- 0.20
    total_w <- weight_ncs + weight_pathway + weight_novel
  }
  if (abs(total_w - 1) > 1e-6)
    stop(sprintf("Weights must sum to 1 (current sum: %.3f)", total_w))

  # ---------- 1. NCS score ----------
  # Use absolute value of NCS for NCS < 0 compounds only.
  # Higher absolute NCS = stronger reversal of disease signature = better candidate.
  # NCS > 0 compounds receive score_ncs = 0.
  results$score_ncs <- ifelse(
    !is.na(results$NCS) & results$NCS < 0,
    abs(results$NCS),
    0
  )

  max_ncs <- max(results$score_ncs, na.rm = TRUE)
  if (is.finite(max_ncs) && max_ncs > 0) {
    results$score_ncs <- results$score_ncs / max_ncs
  } else {
    results$score_ncs <- 0
  }

  # ---------- 2. FDR significance score ----------
  # Converts WTCS_FDR to (0 to 1): lower FDR = higher score.
  # Only NCS < 0 compounds receive FDR credit.
  if ("WTCS_FDR" %in% colnames(results)) {
    results$score_fdr <- ifelse(
      is.na(results$WTCS_FDR),
      0,
      pmin(-log10(results$WTCS_FDR + 1e-300) / 10, 1)
    )
  } else {
    warning("WTCS_FDR column not found; FDR score set to 0")
    results$score_fdr <- 0
  }

  # Remove FDR credit from NCS > 0 compounds
  results$score_fdr <- ifelse(
    !is.na(results$NCS) & results$NCS < 0,
    results$score_fdr,
    0
  )

  # ---------- 3. Connectivity score: NCS (75%) + FDR (25%) ----------
  results$score_connectivity <- 0.75 * results$score_ncs +
    0.25 * results$score_fdr

  # ---------- 4. Pathway keyword dictionaries ----------

  # Conservative pool: 4 classical beta-cell reprogramming pathways
  conservative_pathways <- list(
    Wnt        = c("wnt", "gsk-3", "gsk3", "beta-catenin", "frizzled"),
    TGFb       = c("tgf-beta", "tgfb", "alk inhibit", "alk5",
                   "smad", "activin"),
    cAMP_GLP1  = c("adenylyl cyclase", "phosphodiesterase", "pde inhibit",
                   "pka", "forskolin", "epac", "glp-1", "glp1"),
    Epigenetic = c("hdac inhibit", "dnmt inhibit", "bet inhibit",
                   "bromodomain", "histone methyltransferase",
                   "histone deacetylase", "ezh2",
                   "retinoid receptor", "rar agonist", "rxr")
  )

  # Comprehensive pool: 11 validated pathways
  comprehensive_validated <- c(conservative_pathways, list(
    Notch      = c("notch", "gamma secretase", "gamma-secretase"),
    DYRK1A     = c("dyrk1a", "harmine", "dyrk"),
    PI3K_mTOR  = c("pi3k", "mtor", "akt inhibit", "raptor"),
    BMP        = c("bmp signaling", "bone morphogenetic", "bmpr"),
    Senescence = c("senolytic", "bcl-2 inhibit", "navitoclax", "abt-263",
                   "fisetin"),
    ER_stress  = c("er stress", "chaperone", "perk inhibit",
                   "ire1", "atf6", "tudca"),
    JAK_STAT   = c("jak inhibit", "stat inhibit",
                   "ruxolitinib", "tofacitinib")
  ))

  # Novel pathways: derived from intra-islet C0/C6 communication analysis
  # First reported in this study; not available in any existing tool
  novel_pathways_default <- list(
    MDK_axis     = c("midkine", "syndecan", "nucleolin"),
    Somatostatin = c("somatostatin", "sstr", "octreotide"),
    Opioid_C6C0  = c("opioid", "delta receptor", "oprd",
                     "enkephalin", "pomc"),
    GLP1_intra   = c("glucagon receptor", "glp-1 receptor agonist",
                     "exendin", "liraglutide", "semaglutide")
  )

  # Resolve which dictionaries to use
  validated_pathways <- if (!is.null(custom_validated_pathways)) {
    custom_validated_pathways
  } else if (mode == "conservative") {
    conservative_pathways
  } else {
    comprehensive_validated
  }
  novel_pathways <- if (!is.null(custom_novel_pathways)) {
    custom_novel_pathways
  } else {
    novel_pathways_default
  }

  # ---------- 5. Pathway scoring helper ----------
  # Score: 0 (no hit), 0.5 (1 pathway), 1.0 (2+ pathways)
  score_pathway_fn <- function(moa_str, pathway_list) {
    if (is.na(moa_str) || moa_str == "")
      return(list(score = 0, matched = NA_character_))
    moa_lower    <- tolower(moa_str)
    hit_pathways <- names(Filter(function(kws) {
      any(sapply(kws, function(kw) grepl(kw, moa_lower, fixed = FALSE)))
    }, pathway_list))
    n_hits  <- length(hit_pathways)
    score   <- min(n_hits / 2, 1)
    matched <- if (n_hits > 0) paste(hit_pathways, collapse = "|")
    else NA_character_
    list(score = score, matched = matched)
  }

  # ---------- 6. Apply pathway scoring ----------
  moa_col <- intersect(c("moa", "MOA", "moa_x", "moa.x"), colnames(results))
  if (length(moa_col) == 0) {
    warning("No MOA column found; pathway scores set to 0")
    results$score_pathway   <- 0
    results$matched_pathway <- NA_character_
    if (mode == "comprehensive") {
      results$score_novel   <- 0
      results$matched_novel <- NA_character_
    }
  } else {
    val_res <- lapply(results[[moa_col[1]]],
                      score_pathway_fn, validated_pathways)
    results$score_pathway   <- sapply(val_res, `[[`, "score")
    results$matched_pathway <- sapply(val_res, `[[`, "matched")

    if (mode == "comprehensive") {
      nov_res <- lapply(results[[moa_col[1]]],
                        score_pathway_fn, novel_pathways)
      results$score_novel   <- sapply(nov_res, `[[`, "score")
      results$matched_novel <- sapply(nov_res, `[[`, "matched")
    }
  }

  # ---------- 7. Final reprogramming score ----------
  if (mode == "conservative") {
    results$reprogramming_score <-
      weight_ncs     * results$score_connectivity +
      weight_pathway * results$score_pathway
  } else {
    results$reprogramming_score <-
      weight_ncs     * results$score_connectivity +
      weight_pathway * results$score_pathway +
      weight_novel   * results$score_novel
  }

  # ---------- 8. Beta-cell toxicity blacklist (fuzzy match) ----------
  if (is.null(toxicity_blacklist)) {
    toxicity_blacklist <- character(0)
  }

  default_toxicity <- c(
    "streptozotocin",  # selective beta-cell toxin (T1D model)
    "alloxan",         # selective beta-cell toxin
    "cycloheximide",   # protein synthesis inhibitor
    "doxorubicin",     # anthracycline chemotherapy
    "camptothecin",    # topoisomerase poison
    "colchicine",      # microtubule poison
    "thapsigargin",    # ER stress inducer
    "brefeldin"        # Golgi disruptor
  )
  all_toxicity <- union(default_toxicity, tolower(toxicity_blacklist))

  drug_col <- intersect(c("pert_iname", "pert", "name"), colnames(results))
  if (length(drug_col) > 0) {
    drug_names       <- tolower(results[[drug_col[1]]])
    results$is_toxic <- sapply(drug_names, function(x) {
      any(sapply(all_toxicity, function(tox) grepl(tox, x, fixed = TRUE)))
    })
  } else {
    results$is_toxic <- FALSE
  }

  # ---------- 9. Sort and report ----------
  results$scoring_mode <- mode
  results <- results[order(results$reprogramming_score, decreasing = TRUE), ]

  if (mode == "conservative") {
    message(sprintf(
      "[Conservative] %d compounds scored | Pathway hits: %d | Toxic flagged: %d",
      nrow(results),
      sum(results$score_pathway > 0, na.rm = TRUE),
      sum(results$is_toxic, na.rm = TRUE)
    ))
  } else {
    message(sprintf(
      "[Comprehensive] %d compounds scored | Validated: %d | Novel: %d | Toxic: %d",
      nrow(results),
      sum(results$score_pathway > 0, na.rm = TRUE),
      sum(results$score_novel   > 0, na.rm = TRUE),
      sum(results$is_toxic, na.rm = TRUE)
    ))
  }

  return(results)
}


#' Compare conservative and comprehensive scoring modes
#'
#' Computes Spearman rank correlation between two scoring modes to assess
#' robustness of drug prioritization. High correlation indicates that top
#' candidates are stable regardless of pathway assumptions.
#'
#' @param results_cons Dataframe from score_reprogramming(mode="conservative")
#' @param results_comp Dataframe from score_reprogramming(mode="comprehensive")
#' @param top_n Integer. Number of top candidates to compare (default 20).
#'
#' @return Invisibly returns a list with: rho (Spearman correlation),
#'   pval, overlap (shared top candidates), and a printed summary.
#' @export
#'
#' @examples
#' \dontrun{
#' compare_modes(results_cons, results_comp, top_n = 20)
#' }
compare_modes <- function(results_cons, results_comp, top_n = 20) {

  for (col in c("pert_iname", "reprogramming_score")) {
    if (!col %in% colnames(results_cons))
      stop(sprintf("results_cons missing column: %s", col))
    if (!col %in% colnames(results_comp))
      stop(sprintf("results_comp missing column: %s", col))
  }

  common_drugs <- intersect(results_cons$pert_iname,
                            results_comp$pert_iname)
  if (length(common_drugs) == 0)
    stop("No common drugs found between the two results")

  cons_scores <- results_cons$reprogramming_score[
    match(common_drugs, results_cons$pert_iname)]
  comp_scores <- results_comp$reprogramming_score[
    match(common_drugs, results_comp$pert_iname)]

  sp <- stats::cor.test(cons_scores, comp_scores, method = "spearman")

  top_cons <- head(results_cons$pert_iname, top_n)
  top_comp <- head(results_comp$pert_iname, top_n)
  overlap  <- intersect(top_cons, top_comp)

  message("=== Mode Comparison (Spearman Rank Correlation) ===")
  message(sprintf("  rho = %.3f,  p-value = %.4f", sp$estimate, sp$p.value))
  message(sprintf("  Top-%d overlap: %d / %d", top_n, length(overlap), top_n))
  message("  Shared top candidates: ",
          paste(overlap, collapse = ", "))

  invisible(list(
    rho     = as.numeric(sp$estimate),
    pval    = sp$p.value,
    overlap = overlap
  ))
}
