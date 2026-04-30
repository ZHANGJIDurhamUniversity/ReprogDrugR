#' Query CMap/LINCS database with a gene signature
#'
#' @param signature A list with upgenes and downgenes (from build_signature)
#' @param db Database to query, either "lincs" or "lincs2" (default "lincs")
#' @param n_top Number of top results to return (default 100)
#' @param tau Logical, whether to compute tau score (default FALSE)
#' @param workers Number of parallel workers for database query (default 1)
#' @param add_annotations Logical, whether to merge LINCS compound
#'   annotations (default TRUE)
#'
#' @return A dataframe of top drug candidates with scores and LINCS
#'   annotations
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom utils data
#' @export
query_cmap <- function(signature,
                       db              = "lincs",
                       n_top           = 100,
                       tau             = FALSE,
                       workers         = 1,
                       add_annotations = TRUE) {

  if (!is.list(signature) ||
      !all(c("upgenes", "downgenes") %in% names(signature))) {
    stop("`signature` must be a list containing `upgenes` and `downgenes`.")
  }

  up   <- as.character(signature$upgenes)
  down <- as.character(signature$downgenes)

  # ---------- Gene symbol to Entrez ID conversion ----------
  .convert_genes <- function(genes) {
    if (length(genes) == 0) return(character(0))
    genes <- unique(as.character(genes)[
      !is.na(as.character(genes)) & as.character(genes) != ""])
    if (length(genes) == 0) return(character(0))
    is_numeric <- grepl("^[0-9]+$", genes)
    result <- genes
    symbols_to_convert <- genes[!is_numeric]
    if (length(symbols_to_convert) > 0) {
      converted <- as.character(AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys      = symbols_to_convert,
        column    = "ENTREZID",
        keytype   = "SYMBOL",
        multiVals = "first"
      ))
      result[!is_numeric] <- converted
    }
    unique(result[!is.na(result) & result != ""])
  }

  up   <- .convert_genes(up)
  down <- .convert_genes(down)

  if (length(up) == 0)
    warning("No valid up-regulated genes after conversion.")
  if (length(down) == 0)
    warning("No valid down-regulated genes after conversion.")
  if (length(up) == 0 && length(down) == 0)
    stop("No valid genes after Entrez ID conversion.")

  # ---------- LINCS query ----------
  qsig <- signatureSearch::qSig(
    query       = list(upset = up, downset = down),
    gess_method = "LINCS",
    refdb       = db
  )

  result <- tryCatch(
    signatureSearch::gess_lincs(
      qSig           = qsig,
      sortby         = "NCS",
      tau            = tau,
      addAnnotations = FALSE,
      workers        = workers
    ),
    error = function(e) {
      stop("CMap/LINCS query failed: ", e$message)
    }
  )

  result_df       <- as.data.frame(signatureSearch::result(result))
  result_df       <- head(result_df, n_top)
  result_df$.rank <- seq_len(nrow(result_df))

  # ---------- Merge LINCS annotations ----------
  if (isTRUE(add_annotations)) {
    result_df <- tryCatch({
      annot_env <- new.env(parent = emptyenv())
      data("lincs_pert_info", package = "ReprogDrugR", envir = annot_env)

      if (!exists("lincs_pert_info", envir = annot_env))
        stop("`lincs_pert_info` could not be loaded.")

      lincs_pert_info <- as.data.frame(
        get("lincs_pert_info", envir = annot_env)
      )

      annot_by_name          <- lincs_pert_info
      annot_by_name$join_key <- annot_by_name$pert_iname
      annot_by_id            <- lincs_pert_info
      annot_by_id$join_key   <- annot_by_id$pert_id

      annot_long <- rbind(annot_by_name, annot_by_id)
      annot_long <- annot_long[
        !is.na(annot_long$join_key) & annot_long$join_key != "", ,
        drop = FALSE]
      annot_long <- annot_long[
        !duplicated(annot_long$join_key), , drop = FALSE]

      overlapping_cols <- intersect(
        colnames(result_df),
        setdiff(colnames(annot_long), "join_key")
      )
      if (length(overlapping_cols) > 0)
        annot_long <- annot_long[
          , !colnames(annot_long) %in% overlapping_cols,
          drop = FALSE]

      merged_df <- merge(result_df, annot_long,
                         by.x  = "pert",
                         by.y  = "join_key",
                         all.x = TRUE,
                         sort  = FALSE)
      merged_df <- merged_df[order(merged_df$.rank), , drop = FALSE]
      merged_df$.rank <- NULL
      rownames(merged_df) <- NULL
      merged_df
    },
    error = function(e) {
      warning("LINCS annotation failed: ", e$message)
      result_df$.rank <- NULL
      result_df
    })
  } else {
    result_df$.rank <- NULL
  }

  return(result_df)
}
