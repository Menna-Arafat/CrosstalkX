#' Prepare Expression Matrix for Crosstalk Analysis
#'
#' Subsets and preprocesses a gene expression matrix for downstream crosstalk analysis.
#' The function optionally performs group-aware missingness filtering and imputation,
#' gene identifier conversion, and per-gene z-score normalization.
#'
#' @param expr A numeric matrix or data.frame with **genes in rows** and
#'   **samples in columns**. Gene identifiers must be stored in `rownames(expr)`.
#'
#' @param group_names Optional character vector of group-identifying patterns.
#'   Each pattern is matched against `colnames(expr)` using `grepl()` to assign
#'   samples to groups (e.g., `c("Tumor","Control")`). Required for group-aware
#'   missing value handling.
#'
#' @param genes Character vector of gene identifiers to retain. Only rows whose
#'   row names intersect with `genes` are kept. Order follows the intersection.
#'
#' @param z_norm Logical; if `TRUE` (default), performs per-gene z-score
#'   normalization across samples (row-wise scaling).
#'
#' @param missing_thr Numeric threshold in `[0,1]` controlling filtering of
#'   genes based on missingness (default = 0.6). A gene is removed if:
#'   \itemize{
#'     \item The difference between maximum and minimum group-wise NA proportions
#'           is greater than or equal to `missing_thr`, or
#'     \item The overall NA proportion across all samples is greater than or equal
#'           to `missing_thr`.
#'   }
#'
#' @param convert.ID Optional flag. If not `NULL`, gene identifiers are converted
#'   using `map_me(rownames(expr), species)`. The value of `convert.ID` itself is
#'   not passed to `map_me()`; it serves only as a trigger.
#'
#' @param species Character string passed to `map_me()` for identifier mapping. species supported are "human" and "mmusculus"
#'   Default is `"human"`.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Coerces `expr` to a numeric matrix.
#'   \item Subsets rows to genes present in both `expr` and `genes`.
#'   \item If missing values are present:
#'     \itemize{
#'       \item If `group_names` is `NULL`, a message is displayed and no imputation
#'             is performed.
#'       \item If `group_names` is provided:
#'         \enumerate{
#'           \item Samples are assigned to groups by pattern matching.
#'           \item The function stops if any column fails to match a group.
#'           \item Genes with excessive or imbalanced missingness (based on
#'                 `missing_thr`) are removed.
#'           \item Remaining missing values are imputed within each group using
#'                 the row-wise median across that group's samples. If all values
#'                 in a group are missing for a gene, the imputed value defaults
#'                 to 0.
#'         }
#'     }
#'   \item If `convert.ID` is not `NULL`, row names are mapped using `map_me()`.
#'         Unmapped or empty identifiers are removed. If duplicate identifiers
#'         arise, the gene with the highest mean expression is retained.
#'   \item If `z_norm = TRUE`, per-gene z-scores are computed across samples.
#'         Genes with zero variance yield `NA`, which are replaced with 0.
#'         Values are rounded to two decimals.
#' }
#'
#' @return A numeric matrix with genes in rows and samples in columns,
#'   optionally ID-mapped, missingness-filtered, imputed, and/or z-score normalized.
#'
#' @seealso \code{\link{map_me}}
#'
#' @importFrom matrixStats rowMedians
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#'
#' expr <- data.frame(
#'   Tumor_1   = rnorm(5),
#'   Tumor_2   = rnorm(5),
#'   Control_1 = rnorm(5),
#'   Control_2 = rnorm(5)
#' )
#' rownames(expr) <- paste0("Gene", 1:5)
#' expr[2, 2] <- NA
#'
#' result <- prepare_data(
#'   expr = expr,
#'   group_names = c("Tumor", "Control"),
#'   genes = rownames(expr),
#'   z_norm = TRUE,
#'   missing_thr = 0.6
#' )
#'
#' head(result)
#' }


prepare_data = function(expr = NULL, group_names = NULL, genes= NULL, z_norm= TRUE, missing_thr= .6,
                              convert.ID = NULL, species = "human") {

        #species = match.arg(species)
        expr= as.matrix(expr)
        storage.mode(expr) <- "numeric"
        expr <- expr[intersect(genes, rownames(expr)), , drop = FALSE]

        # Missingness handling (requires group_names)
        if (anyNA(expr) && is.null(group_names)) {
          message('Missing values detected. Please provide group_names that appear in column names.\n',
                  'Example: sample_1_control / sample_2_cancer -> group_names = c("control","cancer")')
        }

        if (anyNA(expr) && !is.null(group_names)) {

          group_dist <- lapply(group_names, \(g) which(grepl(g, colnames(expr))))
          # stop if any group has 0 matched columns
          covered_cols <- unique(unlist(group_dist))
          if (length(covered_cols) < ncol(expr)) {
            missing_cols <- setdiff(seq_len(ncol(expr)), covered_cols)
            stop("Unmatched columns: ", paste(colnames(expr)[missing_cols], collapse = ", "))
          }

          # Filtration
          # Filter out genes with imbalanced missing across group and those missing in more than 60% of samples
          # proportion of missing values in each group
          miss_mat <- as.matrix(sapply(group_dist, \(idx) rowMeans(is.na(expr[, idx, drop = FALSE]))))
          # Address imbalance across groups
          miss_dif <- apply(miss_mat, 1, \(v) max(v) - min(v))
          # Overall NA should be less than 60%, missing difference also less than 60%
          keep <- miss_dif < missing_thr & rowMeans(is.na(expr)) < missing_thr
          expr <- expr[keep, , drop = FALSE]

          # Impute with median within each group
          for (gcols in group_dist) {

              na_rows <- which(rowSums(is.na(expr[, gcols, drop = FALSE])) > 0)
              X <- as.matrix(expr[na_rows, gcols, drop = FALSE])
              na_idx <- which(is.na(X), arr.ind = TRUE)
              if (is.null(dim(na_idx))) na_idx <- matrix(na_idx, ncol = 2)
              if (nrow(na_idx) == 0) next

              med <- matrixStats::rowMedians(X, na.rm = TRUE)
              med[is.na(med)] <- 0
              #map each NA cell to its row’s median
              X[na_idx] <- med[na_idx[, 1]]
              expr[na_rows, gcols] <- X
          }
        }

        # Optional ID conversion
        if (!is.null(convert.ID)) {
          ID= map_me(row.names(expr), species)
          ID[is.na(ID)]= ""
          row.names(expr)= ID
          expr = expr[row.names(expr) != "", , drop = FALSE]
          # Drop duplicated genes
          if (anyDuplicated(rownames(expr)) > 0) {
            m <- rowMeans(expr, na.rm = TRUE)
            ord <- order(m, decreasing = TRUE)
            expr <- expr[ord, , drop = FALSE]
            expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
          }
        }

        # Z-score normalization
        if(isTRUE(z_norm)){
          expr = t(scale(t(expr),center = TRUE, scale = TRUE))
          expr[is.na(expr)] <- 0
          expr= round(expr, 2)
          }
        return(expr)
    }


