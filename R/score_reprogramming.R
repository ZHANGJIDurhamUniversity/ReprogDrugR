#' Score drug candidates for beta-cell reprogramming potential
#'
#' Domain-specific scoring system for beta-cell reprogramming drug screening.
#' Two scoring modes are provided to allow sensitivity analysis:
#' \itemize{
#'   \item \code{"conservative"}: NCS (40\%) + validated pathways (60\%),
#'         using 4 classical beta-cell reprogramming pathways
#'         (Wnt, TGFb, cAMP/GLP-1, epigenetic).
#'   \item \code{"comprehensive"}: NCS (30\%) + validated pathways (40\%) +
#'         novel pathways (30\%), where novel pathways are derived from
#'         intra-islet cell-cell communication analysis (C0 vs C6 subtype
#'         communication breakdown), unique to ReprogDrugR.
#' }
#'
#' @param forward_results Dataframe from query_cmap(); must contain an
#'   'NCS' column and an MOA column (one of 'moa', 'MOA', 'moa_x', 'moa.x').
#' @param mode Scoring mode: "conservative" or "comprehensive".
#' @param weight_ncs Numeric. Weight for the CMap NCS score. Defaults:
#'   0.4 (conservative) or 0.3 (comprehensive).
#' @param weight_pathway Numeric. Weight for the validated-pathway score.
#'   Defaults: 0.6 (conservative) or 0.4 (comprehensive).
#' @param weight_novel Numeric. Weight for the novel-pathway score
#'   (comprehensive mode only). Default: 0.3.
#' @param custom_validated_pathways Optional named list to override the
#'   default validated-pathway keyword dictionary.
#' @param custom_novel_pathways Optional named list to override the default
#'   novel-pathway keyword dictionary.
#' @param toxicity_blacklist Character vector of additional drug names to
#'   flag as toxic to beta cells (merged with internal defaults).
#'
#' @return A dataframe sorted by reprogramming_score (descending) with the
#'   following added columns:
#'   \itemize{
#'     \item \code{score_ncs}: Min-max normalized NCS score [0, 1].
#'     \item \code{score_pathway}: Validated pathway score [0, 1].
#'     \item \code{matched_pathway}: Names of matched validated pathways.
#'     \item \code{score_novel}: Novel pathway score (comprehensive only).
#'     \item \code{matched_novel}: Names of matched novel pathways
#'       (comprehensive only).
#'     \item \code{reprogramming_score}: Final weighted score.
#'     \item \code{is_toxic}: Logical, TRUE if drug is on the beta-cell
#'       toxicity blacklist.
#'     \item \code{scoring_mode}: Records which mode was used.
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Conservative mode (main results, recommended for primary findings)
#' results_cons <- score_reprogramming(results_forward, mode = "conservative")
#'
#' # Comprehensive mode (sensitivity analysis, leverages novel pathways)
#' results_comp <- score_reprogramming(results_forward, mode = "comprehensive")
#'
#' # Compare top10 between modes — overlap indicates robust candidates
#' intersect(head(results_cons$pert_iname, 10),
#'           head(results_comp$pert_iname, 10))
#' }
score_reprogramming <- function(forward_results,
                                mode = c("conservative", "comprehensive"),
                                weight_ncs     = NULL,
                                weight_pathway = NULL,
                                weight_novel   = NULL,
                                custom_validated_pathways = NULL,
                                custom_novel_pathways     = NULL,
                                toxicity_blacklist        = NULL) {

  mode <- match.arg(mode)

  # ---------- Input validation ----------
  if (!is.data.frame(forward_results) || nrow(forward_results) == 0)
    stop("forward_results must be a non-empty dataframe from query_cmap()")
  if (!("NCS" %in% colnames(forward_results)))
    stop("forward_results must contain an 'NCS' column")

  # ---------- Default weights ----------
  # Conservative: NCS 40% + validated pathways 60%
  # Comprehensive: NCS 30% + validated 40% + novel 30%
  if (mode == "conservative") {
    if (is.null(weight_ncs))     weight_ncs     <- 0.40
    if (is.null(weight_pathway)) weight_pathway <- 0.60
    weight_novel <- 0
    total_w <- weight_ncs + weight_pathway
  } else {
    if (is.null(weight_ncs))     weight_ncs     <- 0.30
    if (is.null(weight_pathway)) weight_pathway <- 0.40
    if (is.null(weight_novel))   weight_novel   <- 0.30
    total_w <- weight_ncs + weight_pathway + weight_novel
  }
  if (abs(total_w - 1) > 1e-6)
    stop(sprintf("Weights must sum to 1 (current sum: %.3f)", total_w))

  results <- forward_results

  # ---------- 1. NCS score: min-max normalize to [0, 1] ----------
  ncs_min <- min(results$NCS, na.rm = TRUE)
  ncs_max <- max(results$NCS, na.rm = TRUE)
  ncs_range <- ncs_max - ncs_min
  if (ncs_range == 0) {
    # All NCS values identical — assign neutral 0.5 to avoid NaN
    results$score_ncs <- 0.5
  } else {
    results$score_ncs <- (results$NCS - ncs_min) / ncs_range
  }

  # ---------- 2. Pathway keyword dictionaries ----------
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
  # (4 conservative + 7 supporting validated pathways)
  comprehensive_validated <- c(conservative_pathways, list(
    Notch       = c("notch", "gamma secretase", "gamma-secretase"),
    DYRK1A      = c("dyrk1a", "harmine", "dyrk"),
    PI3K_mTOR   = c("pi3k", "mtor", "akt inhibit", "raptor"),
    BMP         = c("bmp signaling", "bone morphogenetic", "bmpr"),
    Senescence  = c("senolytic", "bcl-2 inhibit", "navitoclax",
                    "abt-263", "fisetin"),
    ER_stress   = c("er stress", "chaperone", "perk inhibit",
                    "ire1", "atf6", "tudca"),
    JAK_STAT    = c("jak inhibit", "stat inhibit",
                    "ruxolitinib", "tofacitinib")
  ))

  # Novel pathways: derived from intra-islet C0/C6 communication analysis
  # — unique to this study, not present in any existing drug screening tool
  novel_pathways_default <- list(
    MDK_axis      = c("midkine", "syndecan", "nucleolin"),
    Somatostatin  = c("somatostatin", "sstr", "octreotide"),
    Opioid_C6C0   = c("opioid", "delta receptor", "oprd",
                      "enkephalin", "pomc"),
    GLP1_intra    = c("glucagon receptor", "glp-1 receptor agonist",
                      "exendin", "liraglutide", "semaglutide")
  )

  # Resolve which pathway dictionaries to use
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

  # ---------- 3. Pathway scoring helper ----------
  # For each drug's MOA string, count how many pathways it hits.
  # Score: 0 (no hit), 0.5 (1 pathway), 1.0 (2+ pathways).
  # This rewards multi-pathway drugs (more biologically plausible).
  score_pathway_fn <- function(moa_str, pathway_list) {
    if (is.na(moa_str) || moa_str == "")
      return(list(score = 0, matched = NA_character_))
    moa_lower <- tolower(moa_str)
    hit_pathways <- names(Filter(function(kws) {
      any(sapply(kws, function(kw) grepl(kw, moa_lower, fixed = FALSE)))
    }, pathway_list))
    n_hits  <- length(hit_pathways)
    score   <- min(n_hits / 2, 1)
    matched <- if (n_hits > 0) paste(hit_pathways, collapse = "|")
    else NA_character_
    list(score = score, matched = matched)
  }

  # ---------- 4. Apply pathway scoring ----------
  # Robust MOA column detection (different LINCS releases use different names)
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
    # Validated pathway scoring
    val_results <- lapply(results[[moa_col[1]]],
                          score_pathway_fn, validated_pathways)
    results$score_pathway   <- sapply(val_results, `[[`, "score")
    results$matched_pathway <- sapply(val_results, `[[`, "matched")

    # Novel pathway scoring (comprehensive mode only)
    if (mode == "comprehensive") {
      nov_results <- lapply(results[[moa_col[1]]],
                            score_pathway_fn, novel_pathways)
      results$score_novel   <- sapply(nov_results, `[[`, "score")
      results$matched_novel <- sapply(nov_results, `[[`, "matched")
    }
  }

  # ---------- 5. Compute final reprogramming score ----------
  if (mode == "conservative") {
    results$reprogramming_score <-
      weight_ncs     * results$score_ncs +
      weight_pathway * results$score_pathway
  } else {
    results$reprogramming_score <-
      weight_ncs     * results$score_ncs +
      weight_pathway * results$score_pathway +
      weight_novel   * results$score_novel
  }

  # ---------- 6. Beta-cell toxicity blacklist ----------
  # These compounds are well-known to damage beta cells; they should
  # never be considered as reprogramming candidates regardless of NCS.
  default_toxicity <- c(
    "streptozotocin",  # selective beta-cell toxin (used to model T1D)
    "alloxan",         # selective beta-cell toxin
    "cycloheximide",   # general protein synthesis inhibitor
    "doxorubicin",     # anthracycline chemotherapy
    "camptothecin",    # topoisomerase poison
    "colchicine",      # microtubule poison
    "thapsigargin",    # ER stress inducer (beta cells especially sensitive)
    "brefeldin"        # Golgi disruptor
  )
  all_toxicity <- union(default_toxicity, tolower(toxicity_blacklist))

  drug_col <- intersect(c("pert_iname", "pert", "name"), colnames(results))
  if (length(drug_col) > 0) {
    results$is_toxic <- tolower(results[[drug_col[1]]]) %in% all_toxicity
  } else {
    results$is_toxic <- FALSE
  }

  # ---------- 7. Sort and report ----------
  results$scoring_mode <- mode
  results <- results[order(results$reprogramming_score,
                           decreasing = TRUE), ]

  if (mode == "conservative") {
    message(sprintf(
      "[Conservative mode] Scored %d compounds | Pathway hits: %d | Toxic flagged: %d",
      nrow(results),
      sum(results$score_pathway > 0, na.rm = TRUE),
      sum(results$is_toxic, na.rm = TRUE)
    ))
  } else {
    message(sprintf(
      "[Comprehensive mode] Scored %d compounds | Validated hits: %d | Novel hits: %d | Toxic flagged: %d",
      nrow(results),
      sum(results$score_pathway > 0, na.rm = TRUE),
      sum(results$score_novel   > 0, na.rm = TRUE),
      sum(results$is_toxic, na.rm = TRUE)
    ))
  }

  return(results)
}
