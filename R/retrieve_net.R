#' Prepare a Gene Interaction Network
#'
#' Prepares a gene–gene interaction edge list for downstream crosstalk analysis.
#' If no network is provided, a default interaction network bundled with
#' \pkg{CrosstalkX} is loaded.
#'
#' @param net A data frame representing an interaction edge list. The first two
#'   columns are assumed to correspond to source and target nodes (e.g., genes or TFs).
#'   Additional columns (e.g., mode of regulation) are allowed. If `NULL`, the
#'   default network dataset `net` from \pkg{CrosstalkX} is used.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Loads the default interaction network from \pkg{CrosstalkX} if `net` is `NULL`.
#'   \item Renames the first two columns to `source` and `target`.
#'   \item Removes edges with missing or empty source/target entries.
#'   \item Removes duplicated source–target pairs while preserving additional columns.
#' }
#'
#' @return A `data.frame` containing a cleaned interaction network with at least
#'   the following columns:
#' \describe{
#'   \item{source}{Source node (gene or transcription factor).}
#'   \item{target}{Target node (gene).}
#' }
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Use the default CrosstalkX interaction network
#' net1 <- retrieve_net()
#'
#' # Provide a custom interaction network
#' net2 <- retrieve_net(custom_net)
#' }


retrieve_net = function(net = NULL) {

  # Load base PPI if not provided
  if (is.null(net)) {
    data("net", package = "CrosstalkX", envir = environment())
    net <- get("net", envir = environment())
    }

  colnames(net)[1:2] = c("source", "target")
  # Clean missing/empty entries
  net = net |> dplyr::filter(complete.cases(net) , source != "", target != "")
  net= dplyr::distinct(net, source, target, .keep_all = TRUE)

  return(net)
}
