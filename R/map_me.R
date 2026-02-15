#' Map identifiers to gene symbols (and symbols to Ensembl gene IDs)
#'
#' @description
#' `map_me()` maps a vector of identifiers to gene symbols (or maps gene symbols to Ensembl gene IDs)
#' using Ensembl resources. Depending on `species` and `convert.ID`, mapping is performed via
#' \pkg{biomaRt} (Ensembl BioMart). For human Ensembl gene IDs, an \pkg{EnsDb} backend may be used
#' if implemented in your code path.
#'
#' @param ID Character vector of identifiers to map (e.g., Ensembl gene IDs, UniProt IDs, Entrez IDs, or gene symbols).
#' @param species Character scalar indicating the organism. Supported values are:
#' \itemize{
#'   \item `"human"` (Ensembl dataset: `hsapiens_gene_ensembl`)
#'   \item `"mmusculus"` (Ensembl dataset: `mmusculus_gene_ensembl`)
#' }
#' @param convert.ID Character scalar specifying the conversion to perform. Supported values are:
#' \itemize{
#'   \item `"ensembl.gene_to_symbol"`: Ensembl gene ID \eqn{\rightarrow} gene symbol
#'   \item `"ensembl.protein_to_symbol"`: Ensembl peptide/protein ID \eqn{\rightarrow} gene symbol
#'   \item `"uniprot_to_symbol"`: UniProt GN ID \eqn{\rightarrow} gene symbol (isoform suffixes like `-1` are removed)
#'   \item `"uniprot_swiss_to_symbol"`: UniProt Swiss-Prot accession \eqn{\rightarrow} gene symbol
#'   \item `"entrezid_to_symbol"`: Entrez gene ID \eqn{\rightarrow} gene symbol
#'   \item `"symbol_to_ensembl"`: gene symbol \eqn{\rightarrow} Ensembl gene ID
#' }
#'
#' @details
#' Mapping is performed using Ensembl BioMart via \code{biomaRt::useMart()} and \code{biomaRt::getBM()}.
#' The returned vector is aligned to the input `ID` (same length and order).
#'
#' For `"uniprot_to_symbol"`, isoform suffixes in the input (e.g., `P12345-1`) are removed prior to querying.
#'
#' For `species = "human"`, gene symbols correspond to HGNC symbols; for `species = "mmusculus"`,
#' symbols correspond to MGI symbols.
#'
#' \strong{Important:} BioMart queries require an internet connection and can fail due to Ensembl downtime
#' or network restrictions. Consider caching results for reproducibility.
#'
#' @return A character vector of mapped identifiers (same length as `ID`). Unmapped entries are `NA`.
#'
#' @seealso \code{\link[biomaRt]{useMart}}, \code{\link[biomaRt]{getBM}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Ensembl gene IDs -> symbols
#' ids <- c("ENSG00000141510", "ENSG00000155657")
#' map_me(ID = ids, convert.ID = "ensembl.gene_to_symbol", species = "human")
#'
#' # UniProt Swiss-Prot -> symbols
#' uniprot_ids <- c("P04637", "P38398")
#' map_me(ID = uniprot_ids, convert.ID = "uniprot_swiss_to_symbol", species = "human")
#'
#' # Symbols -> Ensembl gene IDs
#' map_me(ID = c("TP53", "EGFR"), convert.ID = "symbol_to_ensembl", species = "human")
#' }


map_me <- function( ID = NULL, species ="human" , convert.ID=  NULL) {

          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager") }
          pkgs <- c("ensembldb", "biomaRt")
          for (p in pkgs) {
            if (!requireNamespace(p, quietly = TRUE)) { BiocManager::install(p)  }
          }

          if(species %in% c("human" ,"mmusculus" )){
            stop( 'species must be either "human" or "mmusculus"')
          }
          ID <- as.character(ID)

          # Map using provided mapping file
          # if (!is.null(mapping_file)) {
          #   mapping_df <- as.data.frame(mapping_file)
          #   if (ncol(mapping_df) < 2) stop("mapping_file must have at least 2 columns: query, mapped.")
          #   colnames(mapping_df)[1:2] <- c("query", "gene")
          #   mapped <- mapping_df$gene[match(ID, mapping_df$query)]
          #   return(mapped)
          # }
           # Map EnsDb.Hsapiens.v79
          # if (species== "human" && convert.ID == "ensembl.gene_to_symbol"){
          #   mapping_df <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ID,
          #                                   keytype = "GENEID", columns = c("SYMBOL","GENEID"))
          #   mapped_ID= mapping_df$SYMBOL[match(ID, mapping_df$GENEID)]
          #   return(mapped_ID)
          # }

          # Map biomart
          if (!is.null(convert.ID)) {
            valid_conversions <- c("ensembl.gene_to_symbol", "ensembl.protein_to_symbol", "uniprot_swiss_to_symbol",
                                   "uniprot_to_symbol", "entrezid_to_symbol", "symbol_to_ensembl")
            if (!convert.ID %in% valid_conversions) {
              stop("Invalid 'convert.ID' argument. Valid options are: ",
                   paste(valid_conversions, collapse = ", "))
            }

            if(species== "mmusculus"){
              ensembl<- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
              symbol=  "mgi_symbol"
            }
            if (species== "human" & convert.ID != "ensembl.gene_to_symbol"){
              ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
              symbol= "hgnc_symbol"
            }

            # Process each column
            if (convert.ID == "uniprot_to_symbol"){
              ID=   gsub("-\\d+$", "", ID)
            }

            # Define biomart attributes and filters
            conversion_params <- list(
              "ensembl.gene_to_symbol" = list(attr = c("ensembl_gene_id", symbol), filter = "ensembl_gene_id"),
              "ensembl.protein_to_symbol" = list(attr = c("ensembl_peptide_id", symbol), filter = "ensembl_peptide_id"),
              "uniprot_to_symbol" = list(attr = c("uniprot_gn_id", symbol), filter = "uniprot_gn_id"),
              "uniprot_swiss_to_symbol"= list(attr= c("uniprotswissprot", symbol), filter= "uniprotswissprot"),
              "entrezid_to_symbol" = list(attr = c("entrezgene_id", symbol), filter = "entrezgene_id"),
              "symbol_to_ensembl"= list(attr=    c(symbol, "ensembl_gene_id"),filter=  "hgnc_symbol")
            )

           # Map Biomart
            params <- conversion_params[[convert.ID]]
            values <- unique(ID)
            mapping_df = biomaRt::getBM(
                                    attributes = params$attr,
                                    filters = params$filter,
                                    values = values,
                                    mart = ensembl)
            mapping_df= as.data.frame(mapping_df)
            names(mapping_df)= c("query", "gene")
            mapped_ID= mapping_df$gene[match(ID, mapping_df$query)]

          }

          return(mapped_ID)
}
