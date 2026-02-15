
#' @importFrom stats complete.cases p.adjust sd
#' @importFrom utils data install.packages stack
#' @importFrom dplyr desc
#' @importFrom rlang .data
#' @import pbapply
NULL


utils::globalVariables(c(
  "p_value", "pvalue", "score", "target", "mor",
  "Degree", "edges", "nodes", "weight",
  "x_val", "y_val", "size", "shape",
  "label_size", "label_length", "new_col",
  "pval_Degree", "pval_MI_norm", "name"
))

