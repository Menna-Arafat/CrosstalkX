#' Compute Mutual Information for Gene Pairs in a Network
#'
#' Computes mutual information (MI) for each gene–gene interaction (edge) in a provided
#' network using the k-nearest neighbor (kNN) MI estimator implemented in
#' \pkg{parmigene}. Only edges whose endpoints are present in the expression matrix
#' are evaluated.
#'
#' @param z.score A numeric matrix or data.frame of expression-like values with
#'   **genes as rows** and **samples as columns**. Row names must correspond to gene
#'   identifiers.
#'
#' @param net A data.frame describing gene–gene interactions. Only the first two
#'   columns are used and are interpreted as source (`From`) and target (`To`) genes.
#'
#' @param k Integer. Number of nearest neighbors used by the kNN MI estimator.
#'   Passed to \code{parmigene::knnmi()}. Default is `3`.
#'
#' @param noise Numeric. Small noise value added internally by the estimator to
#'   break ties. Passed to \code{parmigene::knnmi()}. Default is `1e-12`.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Coerces `z.score` to a numeric matrix.
#'   \item Keeps only the first two columns of `net` and renames them to `From` and `To`.
#'   \item Filters the network to retain edges for which both genes are present in
#'   `rownames(z.score)`.
#'   \item Computes mutual information for each retained edge using
#'   \code{parmigene::knnmi()}.
#'   \item Replaces any `NA` MI values with `0`.
#' }
#'
#' If no edges remain after filtering, the function stops with an error.
#'
#' @return A `data.frame` containing the filtered interaction network with an
#' additional column:
#' \describe{
#'   \item{From}{Source gene identifier.}
#'   \item{To}{Target gene identifier.}
#'   \item{mi_score}{Estimated mutual information for the gene pair.}
#' }
#'
#' @seealso \code{\link[parmigene]{knnmi}}
#'
#' @importFrom parmigene knnmi
#'
#' @export
#'
#' @examples
#' \dontrun{
#' z <- matrix(rnorm(1000), nrow = 10)
#' rownames(z) <- paste0("gene", 1:10)
#' colnames(z) <- paste0("s", 1:100)
#'
#' net <- data.frame(
#'   geneA = c("gene1", "gene2"),
#'   geneB = c("gene3", "gene4")
#' )
#'
#' mi_df <- compute_mi(z, net, k = 3)
#' head(mi_df)
#' }

compute_mi <- function(z.score, net, k = 3, noise = 1e-12) {

      x <- as.matrix(z.score)
      net= net[,1:2]
      names(net)= c("From", "To")
      # keep only valid edges
      net2 <- net |> dplyr::filter(net$From %in% rownames(x) & net$To %in% rownames(x))
      if (!nrow(net2)) stop("No edges matched rownames(z.score).")

      net2$mi_score <- vapply(seq_len(nrow(net2)), function(i) {
        parmigene::knnmi(x[net2$From[i], ], x[net2$To[i], ], k = k, noise = noise)
      }, numeric(1))

      net2$mi_score[is.na(net2$mi_score)] <- 0
      return(net2)
}
