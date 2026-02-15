#' Compute Pathway–Pathway Interaction Statistics from a Gene–Gene MI Network
#'
#' Compute pairwise interaction metrics between pathways/terms using (i) a gene→term
#' membership graph and (ii) a gene–gene interaction graph carrying mutual information
#' (MI) weights. For each pathway pair, the function summarizes the number of connecting
#' gene–gene edges (degree) and the total/normalized MI across those edges.
#'
#' @param p_graph A directed \pkg{igraph} object encoding pathway membership edges.
#'   After conversion with [igraph::as_edgelist()], the first column is interpreted as
#'   genes and the second as terms/pathways (i.e., each edge indicates a gene belongs to
#'   a pathway).
#'
#' @param g_graph An undirected \pkg{igraph} gene–gene interaction graph. Edges are expected
#'   to carry a numeric attribute named `mi_score`, used as MI weights via
#'   [igraph::as_adjacency_matrix()]. The unweighted adjacency is used to
#'   compute interaction degree.
#'
#' @param terms Character vector of term/pathway identifiers used as the column names of the
#'   incidence matrix (must match the term labels in `p_graph`).
#'
#' @param genes Character vector of gene identifiers used as the row names of the incidence
#'   matrix (must match the gene labels in `p_graph` and the vertices in `g_graph`).
#'
#' @details
#' The pathway membership graph is converted to a sparse gene × term incidence matrix \eqn{M}.
#' Two gene–gene adjacency matrices are derived from `g_graph`:
#' \itemize{
#'   \item \eqn{A}: weighted adjacency using edge attribute `mi_score`
#'   \item \eqn{B}: unweighted adjacency (presence/absence)
#' }
#'
#' Pathway–pathway matrices are computed as:
#' \deqn{MI = M^\top A M}
#' \deqn{Degree = M^\top B M}
#' Diagonals are set to zero (no self-interactions). The output is restricted to unique
#' pathway pairs (upper triangle).
#'
#' Normalization uses pathway sizes \eqn{|p|} and \eqn{|q|} (column sums of \eqn{M}):
#' \deqn{mi\_norm(p,q) = mi(p,q) / (|p| + |q|)}
#' and similarly for degree (`deg_norm`).
#'
#' @return A `data.frame` with one row per unique pathway pair and columns:
#' \describe{
#'   \item{Term1}{First pathway/term identifier.}
#'   \item{Term2}{Second pathway/term identifier.}
#'   \item{degree}{Number of gene–gene edges connecting the two terms.}
#'   \item{deg_norm}{Degree normalized by \eqn{|p| + |q|}.}
#'   \item{mi_score}{Sum of `mi_score` across gene–gene edges connecting the two terms.}
#'   \item{mi_norm}{`mi_score` normalized by \eqn{|p| + |q|}.}
#' }
#'
#' @importFrom igraph as_edgelist as_adjacency_matrix
#' @importFrom Matrix sparseMatrix t colSums
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' out <- compute_interactions(p_graph, g_graph, terms, genes)
#' head(out)
#' }


compute_interactions <- function(p_graph, g_graph, terms, genes) {

        # Turn pathway graph to gene X pathway adjacency matrix
        edges <- as.data.frame(igraph::as_edgelist(p_graph), stringsAsFactors = FALSE)
        names(edges) <- c("gene", "term")

        # M = gene X pathway adjacency matrix
        M <- Matrix::sparseMatrix(
          i = match(edges[["gene"]], genes),
          j = match(edges[["term"]], terms),
          x = 1,
          dims = c(length(genes), length(terms)),
          dimnames = list(genes, terms)
        )

        # Gene–gene adjacency matrices
        A <- igraph::as_adjacency_matrix(g_graph, attr = "mi_score", sparse = TRUE)
        A[is.na(A)] <- 0
        diag(A) <- 0

        B <- sign(igraph::as_adjacency_matrix(g_graph, sparse = TRUE))
        diag(B) <- 0

        # All pathway–pathway interactions
        mi_mat  <- Matrix::t(M) %*% A %*% M
        deg_mat <- Matrix::t(M) %*% B %*% M
        diag(mi_mat)  <- 0
        diag(deg_mat) <- 0

        # Normalization
        sizes <- Matrix::colSums(M)
        den <- outer(sizes, sizes, "+")
        mi_norm  <- mi_mat / den
        deg_norm <- deg_mat / den

        mi_long      <- reshape2::melt(as.matrix(mi_mat),
                                       varnames = c("Term1", "Term2"),
                                       value.name = "mi_score")
        mi_norm_long <- reshape2::melt(as.matrix(mi_norm),
                                       varnames = c("Term1", "Term2"),
                                       value.name = "mi_norm")
        deg_long     <- reshape2::melt(as.matrix(deg_mat),
                                       varnames = c("Term1", "Term2"),
                                       value.name = "degree")
        deg_norm_long <- reshape2::melt(as.matrix(deg_norm),
                                        varnames = c("Term1", "Term2"),
                                        value.name = "deg_norm")

        final <- Reduce(
          function(x, y) merge(x, y, by = c("Term1", "Term2")),
          list(deg_long, deg_norm_long, mi_long, mi_norm_long)
        )

        # remove self-interactions
        final <- final[final[["Term1"]] != final[["Term2"]], , drop = FALSE]

        # keep upper triangle (unique pairs)
        final <- final[as.character(final[["Term1"]]) < as.character(final[["Term2"]]), , drop = FALSE]
        final <- final[!duplicated(final[c("Term1", "Term2")]), , drop = FALSE]

        final
}

# old
# cross_edges <- igraph::E(graph)[igraph::V(graph)[p1] %--% igraph::V(graph)[p2]]
# mi_score= sum(unlist( igraph::E(graph)[cross_edges]$mi_score), na.rm = TRUE)/ (n1+n2)
