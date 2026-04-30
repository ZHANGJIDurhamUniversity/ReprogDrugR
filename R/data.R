#' LINCS perturbation information
#'
#' Annotation data for LINCS compounds including mechanism of action (MOA),
#' target, clinical phase, and chemical identifiers. Used internally by
#' query_cmap() to annotate drug screening results.
#'
#' @format A data frame with the following key columns:
#' \describe{
#'   \item{pert_iname}{Compound name}
#'   \item{pert_id}{LINCS perturbation ID (BRD number)}
#'   \item{moa}{Mechanism of action}
#'   \item{target}{Primary molecular target}
#'   \item{pubchem_cid}{PubChem compound ID}
#'   \item{chembl_id}{ChEMBL compound ID}
#'   \item{max_phase}{Maximum clinical trial phase}
#' }
#' @source LINCS L1000 database
"lincs_pert_info"
