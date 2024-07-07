library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
#library(SeuratWrappers)
library(monocle3)
library(magrittr)
library(SeuratWrappers)

set.seed(1234)

#Seurat 
lung <- Load10X_Spatial(data.dir = "../demo_data/Spaceranger/")
lung <- subset(lung, subset = nFeature_Spatial > 0)
lung <- SCTransform(lung, assay = "Spatial", verbose = FALSE) %>%
		 	RunPCA( assay = "SCT") %>% 
			FindNeighbors(dims = 1:30) %>% 
			FindClusters() %>% 
			RunUMAP(dims = 1:30)

# chehck Clsuter result
p1 <- SpatialDimPlot(lung, label = TRUE, label.size = 3)
p2 <- DimPlot(lung, reduction = "umap", label = TRUE)
p1 + p2

# Monocle
cds <- as.cell_data_set(lung)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# pseudotime
# root -> cluster 9
get_earliest_principal_node <- function(cds, time_bin = "9"){
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
p3 <- plot_cells(cds, color_cells_by = "seurat_clusters", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
p4 <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

