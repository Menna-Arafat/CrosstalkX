#' Compute Pairwise Jaccard Similarity Between Terms
#'
#' Compute Jaccard similarity between all term pairs based on gene-set overlap.
#' The function converts a gene→term membership graph into a sparse gene × term
#' incidence matrix and derives pairwise intersections/unions via matrix
#' cross-products.
#'
#' @param p_graph A directed \pkg{igraph} object encoding membership edges. After
#'   conversion with [igraph::as_edgelist()], the first column is interpreted as genes
#'   and the second as terms (i.e., each edge indicates that a gene belongs to a term).
#'
#' @param terms Character vector of term identifiers defining the column order
#'   (and column names) of the incidence matrix.
#'
#' @param genes Character vector of gene identifiers defining the row order
#'   (and row names) of the incidence matrix.
#'
#' @details
#' A sparse incidence matrix \eqn{M} is constructed where \eqn{M[g,t] = 1} if gene \eqn{g}
#' belongs to term \eqn{t}. For each term pair \eqn{(t_1, t_2)}:
#' \itemize{
#'   \item Intersections \eqn{|G_{t_1} \cap G_{t_2}|} are obtained from \eqn{M^\top M}.
#'   \item Term sizes \eqn{|G_t|} are obtained from column sums of \eqn{M}.
#'   \item Unions are computed as \eqn{|G_{t_1} \cup G_{t_2}| = |G_{t_1}| + |G_{t_2}| - |G_{t_1} \cap G_{t_2}|}.
#'   \item Jaccard similarity is \eqn{J(t_1,t_2) = |G_{t_1} \cap G_{t_2}| / |G_{t_1} \cup G_{t_2}|}.
#' }
#'
#' The diagonal of the similarity matrix is set to 1, but the returned table excludes
#' self-pairs and retains only unique term pairs (upper triangle).
#'
#' @return A `data.frame` with one row per unique term pair and columns:
#' \describe{
#'   \item{Term1}{First term identifier.}
#'   \item{Term2}{Second term identifier.}
#'   \item{Jaccard}{Jaccard similarity between the two terms.}
#' }
#'
#' @importFrom igraph as_edgelist
#' @importFrom Matrix sparseMatrix t colSums
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#' \dontrun{
#' edges <- data.frame(
#'   gene = c("A","B","C","B","D"),
#'   term = c("T1","T1","T1","T2","T2")
#' )
#' p_graph <- igraph::graph_from_data_frame(edges, directed = TRUE)
#'
#' genes <- unique(edges$gene)
#' terms <- unique(edges$term)
#'
#' J <- compute_jaccard(p_graph, terms, genes)
#' head(J)
#' }


compute_jaccard=function(p_graph, terms, genes) {
        # Turn pathway graph to gene X pathway adjacency matrix
        edges= as.data.frame(igraph::as_edgelist(p_graph), stringsAsFactors = FALSE)
        names(edges)= c("gene", "term")

        # M = gene X pathway adjacency matrix
        M= Matrix::sparseMatrix(i = match(edges$gene, genes),
                                  j = match(edges$term, terms),
                                  x = 1,
                                  dims = c(length(genes), length(terms)),
                                  dimnames = list(genes, terms))

        # Jaccard
        # |Gp ∩ Gq|
        #a matrix of all pairwise intersection
        intersect= Matrix::t(M) %*% M
        sizes= Matrix::colSums(M)
        # |Gp ∪ Gq|
        #a matrix of all pairwise sums
        union= outer(sizes, sizes, "+") - intersect
        J_path= intersect / union
        diag(J_path)= 1

        J= reshape2::melt(as.matrix(J_path), varnames = c("Term1","Term2"), value.name = "Jaccard")
        # remove self-interactions
        J= J[J$Term1 != J$Term2, ]
        # keep upper triangle
        J= J[as.character(J$Term1) < as.character(J$Term2), ]
        J <- J[!duplicated(J[c("Term1","Term2")]), ]

        return(J)
}
