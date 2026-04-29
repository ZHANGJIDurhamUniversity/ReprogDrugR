#' Query CMap/LINCS database in reverse to identify harmful compounds
#'
#' @param signature A list with upgenes and downgenes (from build_signature)
#' @param db Database to query, either "lincs" or "lincs2" (default "lincs")
#' @param n_top Number of top results to return (default 100)
#' @param tau Logical, whether to compute tau score (default FALSE)
#' @param workers Number of parallel workers for database query (default 1)
#'
#' @return A dataframe of compounds that would worsen beta cell state
#' @export
query_cmap_reverse <- function(signature,
                               db = "lincs",
                               n_top = 100,
                               tau = FALSE,
                               workers = 1) {
  # Check signature structure
  if (!is.list(signature) || !all(c("upgenes", "downgenes") %in% names(signature))) {
    stop("`signature` must be a list containing `upgenes` and `downgenes`.")
  }
  # Reverse the signature: swap up and down genes to find harmful compounds
  reversed_signature <- list(
    upgenes = signature$downgenes,
    downgenes = signature$upgenes
  )
  # Reuse query_cmap() with reversed signature
  # All input validation, gene ID conversion and error handling is handled there
  query_cmap(reversed_signature, db = db, n_top = n_top, tau = tau, workers = workers)
}
