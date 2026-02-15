## ----setup, warning=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set( 
  fig.path = "fig/",
  fig.align = "center",
  collapse = TRUE, comment = ">")


## ----warning=FALSE, message=FALSE---------------------------------------------
 # install.packages("devtools")
 # devtools::install_github("Menna-Arafat/CrosstalkX", 
 # build_vignettes = TRUE, dependencies = TRUE)


## ----warning=FALSE, message=FALSE---------------------------------------------
#load lib
library(CrosstalkX)
library(dplyr)
library(tibble)

## ----warning=FALSE, message=FALSE---------------------------------------------
# Expression matrix is expected as genes × samples, with gene IDs in rownames.
# group_names is used downstream to support group-aware preprocessing (e.g., handling NAs).
expr = CrosstalkX::expr |> tibble::column_to_rownames("V1")
group_names = c("Normal", "Tumor")

# Differential expression results are used in two ways:
# 1) fold-change values (linear FC) for pathway enrichment with directionality
# 2) t-statistics for TF activity inference based on downstream targets
de_res = CrosstalkX::de_res |> tibble::column_to_rownames("V1")

# Convert log2 fold change to linear fold change for enrichment inputs that use FC values.
de_res$FC = 2 ^ de_res$log2FoldChange
de_fc = setNames(de_res$FC, row.names(de_res))

# Define DE genes for the gene-prioritization step.
# These represent the "labels" whose local enrichment in the network will be tested.
de_genes = de_res |>
            dplyr::filter(abs(log2FoldChange) >= log2(2), padj <= 0.01) |>
            row.names()

# Compute a signed statistic for TF enrichment.
# The t-statistic captures both effect direction and uncertainty (via the SE).
de_res$t_stat = de_res$log2FoldChange / de_res$lfcSE
tstat_named = setNames(de_res$t_stat, rownames(de_res))

## ----warning=FALSE, message=FALSE---------------------------------------------
# retrieve_net() provides the interaction scaffold used to define each gene's neighborhood.
# If net is not supplied, CrosstalkX loads its bundled default network.
net = retrieve_net()

# prioritize_de() tests each gene of interest for an unusually high number of DE neighbors.
# It uses permutation to build an empirical null distribution, then returns significant genes.
myres = prioritize_de(genes_oi = row.names(expr), de_genes = de_genes, net = net)
top_de= myres$sig
head(top_de, 100)
head(top_de, 100)


## ----warning=FALSE, message=FALSE---------------------------------------------
# Enrichment is performed on prioritized genes rather than on all expressed genes, which
# tends to sharpen signal by focusing on network-coherent regions enriched for DE effects.
gset = CrosstalkX::gset

# Restrict FC vector to prioritized genes before pathway enrichment.
fc_named = de_fc[names(de_fc) %in% row.names(top_de)]

# enrich_path() applies internally either ORA or GSEA enrichment functions from clusterprofiler package.
# it returns enriched pathways with member genes and direction.
# The output is standardized to the columns CrosstalkX expects downstream.
enriched_path = enrich_path(fc_named, enrich = "ORA", gset)
enriched_path = enriched_path |> dplyr::select(ID, geneID, direction)

# TF enrichment uses only regulatory edges (mor != 0) from the network.
# TF activity is inferred from the statistics of downstream targets using functions from decoupleR package.
TF = net |> dplyr::filter(mor != 0)
named_tstat = tstat_named[names(tstat_named) %in% row.names(top_de)]
enriched_tf = enrich_TF(named_tstat, net = TF)

# Combine TF and pathway enrichment into a single term table for crosstalk inference.
# This allows TF–pathway and pathway–pathway couplings to be handled uniformly.
names(enriched_tf) = names(enriched_path)
enrichment_df = rbind(enriched_tf, enriched_path)

# This example focuses on upregulated programs; adjust this filter if you want "down" or both.
enrichment_df = enrichment_df |> dplyr::filter(direction == "up")
row.names(enrichment_df) = NULL
head(enrichment_df)


## ----warning=FALSE, message=FALSE---------------------------------------------
# crosstalk() is the single-call pipeline that:
# prepares term gene sets, preprocesses expression, loads/uses a background network,
# computes MI on gene–gene edges, aggregates MI to term–term scores, and runs permutation tests.
res = crosstalk(enrichment_df, expr, group_names, iter = 100)

## ----warning=FALSE, message=FALSE---------------------------------------------
# filtration:
# low Jaccard overlap reduces trivial links driven by shared genes,
# p-values provide permutation-based evidence of non-random coupling,
# and an MI strength filter helps focus on the stronger half of interactions.
sig = res |>
      dplyr::filter((pval_deg_norm <= 0.05 | pval_mi_norm <= 0.05) & Jaccard < 0.2) |>
      dplyr::filter(mi_norm > quantile(mi_norm, 0.5, na.rm = TRUE))

# A compact string representation of significant couplings for quick inspection.
paste0(gsub("%.*", "", sig$Term1), " >> ", gsub("%.*", "", sig$Term2))


## ----warning=FALSE, message=FALSE---------------------------------------------
# Nodes are labeled with shortened names for readability; colors reflect enrichment direction.
all_labels = unique(c(sig$Term1, sig$Term2))

new_labels = all_labels |>
              stringr::str_sub(1, 50) |>
              gsub("\\%.*|\\*", "", x = _) |>
              gsub("_", " ", x = _)

named_labels = setNames(new_labels, all_labels)

enrichment_df$new_colors = ifelse(enrichment_df$direction == "up", "#CC4248FF", "#2B0B57FF")
idx = match(all_labels, enrichment_df$ID)
new_colors = setNames(enrichment_df$new_colors[idx], all_labels)

# visualize_crosstalk() renders a network view of significant couplings.
visualize_crosstalk(sig, named_labels = named_labels, named_colors = new_colors)


## ----eval=TRUE----------------------------------------------------------------
# scatter_plot() shows scatter plot between jaccard similarity and crosstalk significance. low correlation indicates that inference is not derived by terms similarity..
scatter_plot(sig)


## ----eval=FALSE---------------------------------------------------------------
# 
# # load example data
# enrichment_df= CrosstalkX::enrichment_df
# expr= CrosstalkX::expr |> column_to_rownames("V1")
# 
# # Turn enrichment df to list of pathways. each term has gene hits
# gset_list= prepare_sets(enrichment_df, expr)
# genes= unique(unlist(gset_list))
# 
# # Prepare Expression data
# # prepare_data function applies z-score normalization, impute within each group using median approach, ensure missing values equally distributed between groups
# # Note that it is required to assign group to each sample, e.g. (sample1_Tumor), (sample2_Normal) if there is missing values
# # Optionally map gene IDs to gene symbol using biomart
# group_names= c("Normal", "Tumor")
# z_score= prepare_data(expr, group_names, genes)
# 
# # Compute mutual infrmation between each gene pairs
# mi_df= compute_mi(z_score, net)
# 
# # Compute Crosstalk inference between terms
# res= compute_crosstalk(gset_list, mi_df,iter=100)
# head(res)

