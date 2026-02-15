#' Pathway Enrichment Analysis from Differential Expression Statistics
#'
#' Performs pathway enrichment using either gene set enrichment analysis (GSEA) on a ranked
#' differential-expression vector or over-representation analysis (ORA) on up/down subsets.
#' If `gset` is not supplied, a default gene set collection bundled with \pkg{CrosstalkX}
#' is loaded.
#'
#' @param de_vec A **named** numeric vector of differential expression statistics
#'   (e.g., t-statistics, signed scores, logFC). Names must be gene identifiers compatible
#'   with the gene set collection in `gset`. For GSEA, values are used for ranking.
#'
#' @param enrich Character string specifying the enrichment method. Must be one of
#'   `"GSEA"` or `"ORA"`.
#'
#' @param gset Optional gene set collection. If `NULL`, the dataset `gset` bundled with
#'   \pkg{CrosstalkX} is loaded. Expected format is a named list where each element is a
#'   character vector of gene identifiers; list names correspond to pathway/term IDs.
#'
#' @param pvalueCutoff Numeric. Raw p-value threshold applied after enrichment. Default is `0.05`.
#'
#' @param FDR Numeric. Adjusted p-value (FDR; `p.adjust`) threshold applied after enrichment.
#'   Default is `0.05`.
#'
#' @details
#' The function converts `gset` into a two-column TERM2GENE table and then:
#' \itemize{
#'   \item \strong{GSEA}: ranks `de_vec` in decreasing order and runs
#'   \code{clusterProfiler::GSEA()} with \code{pvalueCutoff = 1} and \code{minGSSize = 1}.
#'   A `direction` column is added based on the sign of NES (`"up"` if NES > 0, else `"down"`).
#'   \item \strong{ORA}: defines two input gene lists, `up` and `down`, and runs
#'   \code{clusterProfiler::enricher()} for each. Results are combined and labeled by `direction`.
#' }
#'
#' Results are filtered using both `pvalueCutoff` and `FDR` thresholds (applied to `pvalue`
#' and `p.adjust`, respectively). The function also adds an `input` column containing the
#' input genes (collapsed with `/`) and writes the results to a CSV file in the current
#' working directory.
#'
#' @return A `data.frame` of enrichment results. The returned columns follow
#' `clusterProfiler` output and typically include pathway/term identifiers, enrichment
#' statistics, `pvalue`, `p.adjust`, and an added `direction` column. An `input` column
#' is appended containing the input gene names collapsed into a single string.
#'
#' @seealso \code{\link[clusterProfiler]{GSEA}}, \code{\link[clusterProfiler]{enricher}}
#'
#' @importFrom dplyr bind_rows
#' @importFrom data.table fwrite
#' @importFrom rlang .data
#' @importFrom clusterProfiler enricher GSEA
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # de_vec must be a named numeric vector (e.g., fold change):
#' de_res$FC= 2 ^ de_res$log2FoldChange
#' de_vec= setNames(de_res$FC, row.names(de_res))
#'
#' gsea_res <- enrich_path(de_vec, enrich = "GSEA", pvalueCutoff = 0.05, FDR = 0.05)
#' head(gsea_res)
#'
#' ora_res  <- enrich_path(de_vec, enrich = "ORA", pvalueCutoff = 0.05, FDR = 0.05)
#' head(ora_res)
#' }

enrich_path= function( de_vec, enrich, gset=NULL, pvalueCutoff= .05, FDR= .05){

        if (is.null(gset)) {
          data("gset", package = "CrosstalkX", envir = environment())
          gset <- get("gset", envir = .GlobalEnv)
        }

          # convert list to long formats
          gmt_long= stack(gset)
          gmt_long= gmt_long[,c(2,1)]

        if (enrich == "GSEA") {
          de_vec <- de_vec[order(de_vec, decreasing= TRUE)]
          e <- clusterProfiler::GSEA(de_vec,
                                     TERM2GENE = gmt_long,
                                     pvalueCutoff = 1,
                                     minGSSize = 1)

          enrich_res <- as.data.frame(e@result)
          enrich_res$direction= ifelse(enrich_res$NES >0, "up", "down")

        } else if (enrich == "ORA") {
          i_up= names(de_vec)[de_vec >1]
          i_down= names(de_vec)[de_vec < 1]
          i_list= list(up= i_up, down= i_down)

          e_list <- lapply(i_list, \(i){
            e=  clusterProfiler::enricher(i,
                                          TERM2GENE = gmt_long,
                                          minGSSize = 1,
                                          pvalueCutoff = 1)
            res <- as.data.frame(e@result)
            res
          })

          enrich_res= dplyr::bind_rows(e_list, .id= "direction")
        } else {
          stop("`enrich` must be either 'GSEA' or 'ORA'.")
        }

        enrich_res= enrich_res |>
                    dplyr::filter(round(.data$p.adjust,2) <= FDR, round(.data$pvalue, 2) <= pvalueCutoff)
        enrich_res$input= paste0(names(de_vec), collapse = "/")

        data.table::fwrite(enrich_res, paste0("enrichment_",enrich, "_pval_", pvalueCutoff,"_FDR_", FDR ,".csv"))
        return(enrich_res)
}

