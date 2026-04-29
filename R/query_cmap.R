#' Query CMap/LINCS database with a gene signature
#'
#' @param signature A list with upgenes and downgenes (from build_signature)
#' @param db Database to query, either "lincs" or "lincs2" (default "lincs")
#' @param n_top Number of top results to return (default 100)
#' @param tau Logical, whether to compute tau score (default FALSE)
#'
#' @return A dataframe of top drug candidates with scores
#' @export
query_cmap <- function(signature, db = "lincs", n_top = 100, tau = FALSE) {

  # Check signature structure
  if (!is.list(signature) || !all(c("upgenes", "downgenes") %in% names(signature))) {
    stop("`signature` must be a list containing `upgenes` and `downgenes`.")
  }

  up <- as.character(signature$upgenes)
  down <- as.character(signature$downgenes)

  # Convert gene symbols to Entrez IDs if necessary
  # Numeric strings are assumed to be Entrez IDs and kept as-is
  # Non-numeric strings are treated as HGNC symbols and converted
  .convert_genes <- function(genes) {
    if (length(genes) == 0) return(character(0))
    is_numeric <- grepl("^[0-9]+$", genes)
    result <- genes
    symbols_to_convert <- genes[!is_numeric]
    if (length(symbols_to_convert) > 0) {
      converted <- as.character(AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = symbols_to_convert,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      ))
      result[!is_numeric] <- converted
    }
    # Remove NAs, empty strings, and deduplicate
    result <- unique(result[!is.na(result) & result != ""])
    return(result)
  }

  up <- .convert_genes(up)
  down <- .convert_genes(down)

  # Warn if either gene set is empty
  if (length(up) == 0) {
    warning("No valid up-regulated genes after conversion.")
  }
  if (length(down) == 0) {
    warning("No valid down-regulated genes after conversion.")
  }

  # Stop if both gene sets are empty
  if (length(up) == 0 && length(down) == 0) {
    stop("No valid genes after Entrez ID conversion. Please check gene symbols or Entrez IDs.")
  }

  # Build query signature
  qsig <- signatureSearch::qSig(
    query = list(upset = up, downset = down),
    gess_method = "LINCS",
    refdb = db
  )

  # Run LINCS query
  result <- tryCatch(
    signatureSearch::gess_lincs(qSig = qsig, sortby = "NCS", tau = tau),
    error = function(e) {
      stop(
        "CMap/LINCS query failed. Please check whether the LINCS reference database is installed and accessible. Original error: ",
        e$message
      )
    }
  )

  result_df <- signatureSearch::result(result)
  return(head(result_df, n_top))
}
