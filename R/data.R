#' Background gene interaction network
#'
#' A background gene–gene interaction network used by CrosstalkX. The network integrates
#' transcription factor–target interactions and protein–protein interactions and is used
#' for mutual information estimation and crosstalk inference.
#'
#' @format A data frame with three columns:
#' \describe{
#'   \item{source}{Source gene or transcription factor}
#'   \item{target}{Target gene}
#'   \item{mor}{Mode of regulation (numeric)}
#' }
#'
#' @source
#' Transcription factor–target interactions from the \pkg{decoupleR} package resources and
#' protein–protein interactions downloaded from STRINGdb, filtered to retain intermediate
#' to high-confidence interactions (combined score > 400).
#'
#' @usage data(net)
"net"


#' Pathway gene sets
#'
#' A collection of pathway gene sets used for pathway enrichment and crosstalk analysis.
#'
#' @format A named list (GMT-style) where each element corresponds to a pathway term and
#'   contains a character vector of gene identifiers.
#'
#' @source
#' Gene sets downloaded from the Bader Lab pathway collection.
#'
#' @usage data(gset)
"gset"


#' Example gene expression data
#'
#' An example gene expression dataset provided for demonstrating CrosstalkX workflows.
#'
#' @format A numeric matrix or data frame with genes as rows and samples as columns.
#'
#' @source
#' The Cancer Genome Atlas (TCGA) project.
#'
#' @usage data(expr)
"expr"


#' Example differential expression results
#'
#' Example results of differential expression analysis used to demonstrate gene
#' prioritization functionality in CrosstalkX.
#'
#' @format A data frame containing gene-level differential expression statistics.
#'
#' @usage data(de_res)
"de_res"


#' Example enrichment input table
#'
#' An example enrichment results table used as input for pathway and transcription factor
#' crosstalk analysis.
#'
#' @format A data.frame containing at least two columns. The function
#'   uses the first two columns only:
#'   \describe{
#'     \item{Column 1}{Pathway/term identifiers (names of gene sets).}
#'     \item{Column 2}{Gene hits as a single character string, with genes separated by
#'       commas (`,`) or slashes (`/`).}
#'   }
#'
#' @usage data(enrichment_df)
"enrichment_df"

#' @author Menna Arafat

