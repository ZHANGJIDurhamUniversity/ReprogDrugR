#' Build a subtype-based gene signature
#'
#' @param highfunc_file Path to the high-functioning beta cell markers CSV
#' @param senescent_file Path to the senescent beta cell markers CSV
#' @param gene_col Column name for gene symbols (default "gene")
#'
#' @return A list with upgenes (C0 markers) and downgenes (C6 markers)
#' @export
build_subtype_signature <- function(highfunc_file,
                                    senescent_file,
                                    gene_col = "gene_symbol") {

  # 读取C0高功能亚群标志基因
  highfunc <- read.csv(highfunc_file)
  upgenes <- highfunc[[gene_col]]

  # 读取C6衰老亚群标志基因
  senescent <- read.csv(senescent_file)
  downgenes <- senescent[[gene_col]]

  # 去除空值
  upgenes <- upgenes[!is.na(upgenes) & upgenes != ""]
  downgenes <- downgenes[!is.na(downgenes) & downgenes != ""]

  return(list(
    upgenes = upgenes,
    downgenes = downgenes
  ))
}
