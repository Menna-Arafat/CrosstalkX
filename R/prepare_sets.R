#' Prepare Pathway Gene Sets from Enrichment Results
#'
#' Converts an enrichment results table (terms and gene hits) into a named list of
#' pathway gene sets, optionally filtering each gene set to genes present in an
#' expression matrix.
#'
#' @param enrichment_df A data.frame containing at least two columns. The function
#'   uses the first two columns only:
#'   \describe{
#'     \item{Column 1}{Pathway/term identifiers (names of gene sets).}
#'     \item{Column 2}{Gene hits as a single character string, with genes separated by
#'       commas (`,`) or slashes (`/`).}
#'   }
#'
#' @param expr A numeric matrix or data.frame of expression values with gene identifiers
#'   in `rownames(expr)`. Gene sets are filtered to retain only genes present in
#'   `rownames(expr)`.
#'
#' @param sep Optional character separator used in the gene-hit column. If `NULL`
#'   (default), the separator is inferred: if any `hits` contains `","` then `","`
#'   is used, otherwise `"/"` is used.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Keeps only the first two columns of `enrichment_df` and renames them to
#'   `terms` and `hits`.
#'   \item Determines a separator (`sep`) if not provided.
#'   \item Converts the term/hits table into a named list using \code{tibble::deframe()}.
#'   \item Splits each hit string into a character vector of genes using \code{strsplit()}.
#'   \item Filters genes to those present in \code{rownames(expr)}.
#'   \item Drops pathways with zero remaining genes.
#' }
#'
#' @return A named list where each element is a character vector of genes for one
#'   pathway/term. Pathways with no genes after filtering are removed.
#'
#' @importFrom tibble deframe
#'
#' @export
#'
#' @examples
#' \dontrun{
#' enrichment_df <- data.frame(
#'   pathwayID = c("Pathway_A", "Pathway_B"),
#'   geneID = c("TP53/BRCA1/EGFR", "MTOR/PIK3CA")
#' )
#'
#' expr <- matrix(rnorm(25), nrow = 5)
#' rownames(expr) <- c("TP53", "BRCA1", "EGFR", "MTOR", "PIK3CA")
#'
#' gset_list <- prepare_sets(enrichment_df, expr)
#' names(gset_list)
#' gset_list[[1]]
#' }

prepare_sets= function(enrichment_df, expr, sep=NULL ){

        enrichment_df=enrichment_df[,1:2]
        names(enrichment_df)[1:2]= c("terms", "hits")
        # Determine separator by inspecting first row
        if(is.null(sep)){
        sep= if (any(grepl(",", enrichment_df$hits))) "," else "/"
        }

        gset_list <- enrichment_df |>
                      tibble::deframe() |>
                      lapply(\(x) strsplit(x, sep, fixed = TRUE)[[1]])

        gset_list <- lapply(gset_list, \(l) l[l %in% row.names(expr)])
        gset_list <- gset_list[lengths(gset_list) > 0]

        return(gset_list)
}

