#' Transcription Factor Enrichment Using decoupleR ULM
#'
#' Infers enriched transcription factors (TFs) from a named vector of gene-level statistics
#' using decoupleR's univariate linear model (ULM). Significant TFs are selected by p-value,
#' optionally restricted to positively scored ("active") regulators, and the corresponding
#' TF \eqn{\rightarrow} target network is returned.
#'
#' @param tstat_vec A **named** numeric vector of gene-level statistics (e.g., t-statistics,
#'   signed scores). Names must be gene identifiers matching the `target` column in `net`.
#'
#' @param net Optional regulatory network (edge list) used by decoupleR. Must contain at least
#'   three columns corresponding to regulator (`source`), target gene (`target`), and mode of
#'   regulation (`mor`). If `NULL`, the dataset `net` bundled with \pkg{CrosstalkX} is loaded
#'   and filtered to `mor != 0`.
#'
#' @param active Logical. If `TRUE` (default), retains only TFs with positive ULM score
#'   (`score > 0`) in addition to `p_value <= pvaluecutoff`. If `FALSE`, retains TFs based
#'   on p-value only.
#'
#' @param pvaluecutoff Numeric. P-value threshold applied to decoupleR results (`p_value`).
#'   Default is `0.05`.
#'
#' @param minGSSize Integer. Minimum number of target genes per TF used by decoupleR ULM
#'   (`minsize` argument). Default is `5`.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Loads the default network from \pkg{CrosstalkX} if `net` is `NULL`, and keeps
#'   only signed edges (`mor != 0`).
#'   \item Converts `tstat_vec` to a matrix and runs \code{decoupleR::run_ulm()} using
#'   `source`, `target`, and `mor` columns from `net`.
#'   \item Filters significant TFs by `p_value <= pvaluecutoff` and, if `active = TRUE`,
#'   additionally requires `score > 0`.
#'   \item Builds a TF \eqn{\rightarrow} target summary table for significant TFs and returns it.
#' }
#'
#' The function writes two CSV files to the current working directory:
#' \itemize{
#'   \item \code{decoupleR_res.csv}: full ULM results from decoupleR.
#'   \item \code{net_target_net.csv}: summarized TF-target network for significant TFs.
#' }
#'
#' @return A `data.frame` with one row per significant TF, containing:
#' \describe{
#'   \item{source}{TF identifier.}
#'   \item{all_targets}{Slash-separated targets (`target`) regulated by the TF among input genes.}
#'   \item{direction}{`"up"` if ULM score > 0, otherwise `"down"`.}
#' }
#'
#' @seealso \code{\link[decoupleR]{run_ulm}}
#'
#' @importFrom dplyr group_by summarise
#' @importFrom data.table fwrite
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # tstat_vec must be a named numeric vector of gene statistics:
#' de_res$t_stat <- de_res$log2FoldChange / de_res$lfcSE
#' # tstat_vec <- setNames(de_res$t_stat, rownames(de_res))
#'
#' tf_targets <- enrich_TF(tstat_vec, net = NULL, active = TRUE, pvaluecutoff = 0.05)
#' head(tf_targets)
#' }


enrich_TF= function(tstat_vec, net, active= TRUE,  pvaluecutoff= .05, minGSSize=5){


        # Load base PPI if not provided
        if (is.null(net)) {
          data("net", package = "CrosstalkX", envir = environment())
          net <- get("net", envir = environment())
          net= net |> dplyr::filter(mor != 0)
        }

        de_mat= as.matrix(tstat_vec)
        # net= get_collectri(organism='human', split_complexes=FALSE)
        # Run Univariate linear model (ulm) and get positive coefficient for active net
        contrast_acts = decoupleR::run_ulm(mat=de_mat,
                                           net=net, .source='source', .target='target',
                                           .mor='mor', minsize = minGSSize)

        # Active nets are those with positive coefficients
        if(isTRUE(active)){
          sig_contrast= contrast_acts |> dplyr::filter(p_value <= pvaluecutoff & score > 0)
        }else{
          sig_contrast= contrast_acts |> dplyr::filter(p_value <= pvaluecutoff)
        }
        sig_contrast$direction= ifelse(sig_contrast$score >0, "up", "down")

        # net target network for active net
        net = net |> dplyr::filter(source %in% sig_contrast$source, target %in% names(tstat_vec))
        net_sig= net |>
                  dplyr::group_by(source) |>
                  dplyr::summarise(all_targets = paste(target, collapse = "/"))
        net_sig$direction= sig_contrast$direction[match(net_sig$source, sig_contrast$source)]

        data.table::fwrite(contrast_acts, "decoupleR_res.csv", row.names = F)
        data.table::fwrite(net_sig, "net_target_net.csv", row.names = F)
        return(net_sig)
}
