#' Build a gene signature from DEG results
#'
#' @param deg_file Path to the DEG CSV file
#' @param lfc_col Column name for log2FoldChange
#' @param fdr_col Column name for FDR/padj
#' @param lfc_cutoff Minimum absolute log2FoldChange (default 1)
#' @param fdr_cutoff Maximum FDR (default 0.05)
#' @param top_n Maximum number of genes per direction (default 150)
#'
#' @return A list with upgenes and downgenes
#' @export
build_signature <- function(deg_file,
                            lfc_col = "log2FoldChange",
                            fdr_col = "padj",
                            lfc_cutoff = 1,
                            fdr_cutoff = 0.05,
                            top_n = 150) {

  deg <- read.csv(deg_file)

  deg <- deg[!is.na(deg$gene_symbol) & deg$gene_symbol != "", ]

  # 过滤显著基因
  deg <- deg[!is.na(deg[[fdr_col]]), ]
  deg <- deg[deg[[fdr_col]] < fdr_cutoff, ]
  deg <- deg[abs(deg[[lfc_col]]) > lfc_cutoff, ]

  # 按log2FC排序取top
  up <- deg[deg[[lfc_col]] > 0, ]
  up <- up[order(up[[lfc_col]], decreasing = TRUE), ]
  upgenes <- head(up$gene_symbol, top_n)

  down <- deg[deg[[lfc_col]] < 0, ]
  down <- down[order(down[[lfc_col]], decreasing = FALSE), ]
  downgenes <- head(down$gene_symbol, top_n)

  return(list(
    upgenes = upgenes,
    downgenes = downgenes
  ))
}
