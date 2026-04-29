#' Query CMap/LINCS database with a gene signature
#'
#' @param signature A list with upgenes and downgenes (from build_signature)
#' @param db Database to query, either "lincs" or "lincs2" (default "lincs")
#' @param n_top Number of top results to return (default 100)
#' @param tau Logical, whether to compute tau score (default FALSE)
#' @param workers Number of parallel workers for database query (default 1)
#' @param add_annotations Logical, whether to add LINCS compound annotations (default TRUE)
#'
#' @return A dataframe of top drug candidates with scores and optional LINCS annotations
#' @export
query_cmap <- function(
    signature,
    db = "lincs",
    n_top = 100,
    tau = FALSE,
    workers = 1,
    add_annotations = TRUE
) {

  # Check signature structure
  if (!is.list(signature) || !all(c("upgenes", "downgenes") %in% names(signature))) {
    stop("`signature` must be a list containing `upgenes` and `downgenes`.")
  }

  up <- as.character(signature$upgenes)
  down <- as.character(signature$downgenes)

  # Convert gene symbols to Entrez IDs if necessary
  .convert_genes <- function(genes) {
    if (length(genes) == 0) return(character(0))

    genes <- as.character(genes)
    genes <- unique(genes[!is.na(genes) & genes != ""])

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

    result <- unique(result[!is.na(result) & result != ""])
    return(result)
  }

  up <- .convert_genes(up)
  down <- .convert_genes(down)

  if (length(up) == 0) {
    warning("No valid up-regulated genes after conversion.")
  }

  if (length(down) == 0) {
    warning("No valid down-regulated genes after conversion.")
  }

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
    {
      signatureSearch::gess_lincs(
        qSig = qsig,
        sortby = "NCS",
        tau = tau,
        addAnnotations = FALSE,
        workers = workers
      )
    },
    error = function(e) {
      stop(
        "CMap/LINCS query failed. Please check whether the LINCS reference database is installed and accessible. Original error: ",
        e$message
      )
    }
  )

  result_df <- signatureSearch::result(result)
  result_df <- head(result_df, n_top)
  result_df <- as.data.frame(result_df)

  # Add LINCS compound annotations
  if (isTRUE(add_annotations)) {

    if (!requireNamespace("signatureSearchData", quietly = TRUE)) {
      warning("Package `signatureSearchData` is not installed. Returning unannotated CMap results.")
      return(result_df)
    }

    annot_env <- new.env()

    annotation_status <- tryCatch(
      {
        data("lincs_pert_info", package = "signatureSearchData", envir = annot_env)

        if (!exists("lincs_pert_info", envir = annot_env)) {
          stop("`lincs_pert_info` could not be loaded from signatureSearchData.")
        }

        lincs_pert_info <- annot_env$lincs_pert_info

        if (!is.data.frame(lincs_pert_info)) {
          stop("`lincs_pert_info` is not a data frame.")
        }

        if (!"pert_iname" %in% colnames(lincs_pert_info)) {
          stop("Column `pert_iname` was not found in `lincs_pert_info`.")
        }

        lincs_pert_info <- as.data.frame(lincs_pert_info)

        # Remove duplicated compound names before merging
        lincs_pert_info <- lincs_pert_info[!duplicated(lincs_pert_info$pert_iname), ]

        # Avoid duplicated column names before merge
        overlapping_cols <- intersect(
          colnames(result_df),
          setdiff(colnames(lincs_pert_info), "pert_iname")
        )

        if (length(overlapping_cols) > 0) {
          lincs_pert_info <- lincs_pert_info[
            ,
            !colnames(lincs_pert_info) %in% overlapping_cols,
            drop = FALSE
          ]
        }

        result_df <- merge(
          result_df,
          lincs_pert_info,
          by.x = "pert",
          by.y = "pert_iname",
          all.x = TRUE,
          sort = FALSE
        )

        result_df
      },
      error = function(e) {
        warning(
          "LINCS annotation failed. Returning unannotated CMap results. Original error: ",
          e$message
        )
        result_df
      }
    )

    result_df <- annotation_status
  }

  return(result_df)
}
