#' Visualize term Cross-talk as a Network Graph
#'
#' Create an undirected network visualization of inferred term–term cross-talk.
#' Nodes represent terms (e.g., pathways/TF programs) and edges represent interactions.
#' Node size is scaled by node degree, and edge width is scaled by the normalized
#' mutual-information score (`mi_norm`).
#'
#' @param sig A `data.frame` of term–term interactions (e.g., output of
#'   [compute_crosstalk()] or [crosstalk()]). Must contain at least the columns
#'   `Term1`, `Term2`, and `mi_norm`. Only these three columns are used for graph
#'   construction.
#'
#' @param named_labels Optional named character vector used to relabel node names via
#'   [dplyr::recode()]. Names correspond to original term IDs and values correspond to
#'   the desired display labels.
#'
#' @param named_colors Optional named character vector of colors (e.g., hex strings)
#'   used to color nodes. Names should match original term IDs (before relabeling).
#'   Nodes without a match (or when `named_colors = NULL`) are assigned a default color.
#'
#' @param top Numeric in (0, 1]. Fraction of highest-degree nodes to retain for visualization.
#'   The function keeps the top `ceiling(n_nodes * top)` nodes by degree. Default is `1`
#'   (retain all nodes).
#'
#' @param layout Character. Layout string passed to [ggraph::ggraph()] (e.g., `"circle"`,
#'   `"fr"`, `"kk"`). Default is `"circle"`.
#'
#' @param prefix Character. Prefix used for output filenames written to disk.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Builds an undirected \pkg{igraph} object from `sig[, c("Term1","Term2","mi_norm")]`.
#'   \item Computes node degree on the full graph and writes a degree table to
#'   `"<prefix>_crosstalkx_terms_degree.csv"`.
#'   \item Retains only the top fraction of nodes by degree as specified by `top`.
#'   \item Assigns node colors using `named_colors` if provided; otherwise uses a default color.
#'   \item Optionally relabels nodes using `named_labels`.
#'   \item Converts to a \pkg{tidygraph} object, rescales node sizes and edge widths, and
#'   renders the network using \pkg{ggraph}.
#'   \item Saves a PNG to `"<prefix>_crossTalkX_plot_top<top>.png"`.
#' }
#'
#' @return A \pkg{ggplot2} object representing the network plot. The plot is also saved
#'   to disk as a PNG. A node-degree table is saved as a CSV.
#'
#' @seealso [compute_crosstalk()], [crosstalk()]
#'
#' @importFrom igraph graph_from_data_frame degree V induced_subgraph
#' @importFrom tidygraph as_tbl_graph activate
#' @importFrom dplyr arrange mutate recode
#' @importFrom scales rescale
#' @importFrom data.table fwrite
#' @import ggraph
#' @import ggplot2
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' # sig is a data.frame with Term1, Term2, and mi_norm
#' p <- visualize_crosstalk(sig, top = 1, layout = "circle")
#' p
#'
#' # optional relabeling and coloring
#' labs <- c("hsa04110" = "Cell cycle", "hsa04010" = "MAPK signaling")
#' cols <- c("hsa04110" = "#1b9e77", "hsa04010" = "#d95f02")
#' p2 <- visualize_crosstalk(sig, named_labels = labs, named_colors = cols, top = 0.5)
#' p2
#' }


visualize_crosstalk= function(sig, named_labels=NULL, named_colors= NULL,
                              top= 1, layout = 'circle', prefix= "res"){

          # Construct graph
          graph= igraph::graph_from_data_frame(sig[, c("Term1","Term2","mi_norm")],  directed = FALSE)
          degree_df= data.frame(term= igraph::V(graph)$name,
                                degree= igraph::degree(graph, mode = "all") ) |>
                                dplyr::arrange(desc(.data$degree), .keep =TRUE)

          # plot top terms
          # top in (0,1]; top = 1 keeps all nodes
          all_degree= igraph::degree(graph, mode = "all")
          k <- max(1L, ceiling(length(all_degree) * top))
          keep_nodes <- names(sort(all_degree, decreasing = TRUE))[seq_len(k)]
          graph <- igraph::induced_subgraph(graph, vids = keep_nodes)

          # color based on direction if enrichment_df has direction column
          if(!is.null(named_colors)){
            igraph::V(graph)$new_col= named_colors[igraph::V(graph)$name]
          }else{
            igraph::V(graph)$new_col= "#CC4248FF"
          }

          # relabel if exists
          if (!is.null(named_labels)) {
            #recode(c("A", "B", "C"), A = "Apple", B = "Banana")
            igraph::V(graph)$name <- dplyr::recode(igraph::V(graph)$name, !!!named_labels)
          }
          # new_label= gsub("\\%.*|\\*", "", igraph::V(graph)$name)
          # new_label= gsub("_", " ", new_label)
          #igraph::V(graph)$name <- stringr::str_sub(new_label, 1, 50)
          #Set Colors for Igraph
          #colors = colorRampPalette(c( "#E3B31C",  "#CC4248FF","#AD305DFF", "#2B0B57FF"))(length(igraph::V(graph))) #"#F1711FFF", "#E35933FF","tan",
          #col = scales::col_numeric(palette = colors , domain = NULL)(degree),

           graph_tbl = tidygraph::as_tbl_graph(graph) |>
                        tidygraph::activate(nodes) |>
                        dplyr::mutate(
                          # Set size based on degree
                          degree =igraph::degree(graph, mode = "all"),
                          size = scales::rescale(.data$degree, to = c(2, 10)),
                          label_size = 4,
                          shape =rep("1", length(igraph::V(graph)$name)),
                          label= igraph::V(graph)$name,
                          label_length = nchar(.data$label)
                        )

          #order nodes by label length
          graph_tbl = graph_tbl |>
                      dplyr::arrange(.data$label_length)

          graph_tbl =graph_tbl |>
                      tidygraph::activate(edges) |>
                      dplyr::mutate(
                        weight = as.numeric(.data$mi_norm),
                        width = scales::rescale(.data$weight, to = c(0.2, 3)))

          #plot
          plot = ggraph::ggraph(graph_tbl, layout = layout) +
                    ggplot2::coord_fixed(expand = TRUE) +
                    #ggplot2::geom_edge_bundle_force(color = "grey55", width= .5)+
                    ggraph::geom_edge_fan(ggplot2::aes(width = .data$width, alpha = 0.7), color = "grey55", show.legend = FALSE) +
                    ggraph::geom_node_point(ggplot2::aes(size = .data$size, color = .data$new_col, shape = .data$shape), show.legend = FALSE) +
                    ggplot2::scale_color_identity() +
                    ggplot2::scale_size_identity() +
                    #ggraph::scale_edge_width(range = c(0.2, 3)) +
                    ggraph::geom_node_text(ggplot2::aes(label = .data$label, size = .data$label_size), color = 'black' , fontface = "bold",
                                                   show.legend = FALSE, repel = TRUE, check_overlap = FALSE,   max.overlaps = Inf) +
                    ggraph::theme_graph(fg_text_colour = 'black' ) + #'black' "white"
                    ggplot2::theme( ggplot2::element_text(size = 16),
                                    axis.ticks = ggplot2::element_blank(),
                                    axis.text = ggplot2::element_blank(),
                                    panel.grid = ggplot2::element_blank(),
                                    panel.border = ggplot2::element_blank(),
                                    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
                                    plot.background = ggplot2::element_rect(fill = "white", colour = NA)
                                    )

          #Save
          ggplot2::ggsave(paste0(prefix, "_crossTalkX_plot_top", top, ".png"), plot, width = 10, height = 10, dpi = 600)
          data.table::fwrite(degree_df, paste0(prefix,"_crosstalkx_terms_degree.csv", row.names = F))

          return(plot)
}
