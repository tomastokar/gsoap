#' Pathways Enrichment of the LUAD Paradoxical Genes
#'
#' The data contains results of the pathway over-representation enrichment
#' analysis of the paradoxical genes obtained from pathDIP (Rahmati et al, 2016).
#' These genes were found deregulated in lung adenocarcinoma (LUAD), and were termed
#' paradoxical because their copy number aberrations are in contrast to
#' their differential expression status (Tokar et al., 2018).
#'
#' @docType data
#'
#' @usage data(pxgenes)
#'
#' @format a data frame with 170 rows and 5 columns:
#' \itemize{
#' \item Source -- original source of the pathway
#' \item Pathway -- pathway name
#' \item p.value -- significance of the pathway enrichment
#' \item FDR -- FDR-adjusted p-value
#' \item Members -- list of query genes members of the pathway
#' }
#'
#' @note The results of the pathway enrichment analysis in this data set
#' and in the original paper may differ, due to the different pathDIP
#' version used.
#'
#' @source
#' \itemize{
#' \item Rahmati et al.
#' "pathDIP: an annotated resource for known and predicted human
#' gene-pathway associations and pathway enrichment analysis."
#' Nucleic acids research 45.D1 (2016): D419-D426.
#'
#' \item Tokar et al.
#' "Differentially expressed microRNAs in lung adenocarcinoma
#' invert effects of copy number aberrations of prognostic genes."
#' Oncotarget 9.10 (2018): 9137.
#' }
#'
#' @examples
#' data(pxgenes)
#' summary(pxgenes)
#'
"pxgenes"
