library(Seurat) # v4.0.0 or 4.3.0
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

# Please input the output directory from Space Ranger
lung <- Load10X_Spatial(data.dir = "<INPUT>")

# number of spots
ncol(lung@assays$Spatial@counts)
#lung <- subset(lung, subset = nFeature_Spatial > 0)
#ncol(lung@assays$Spatial@counts)

# QC plot
pdf("Visium_nCount.pdf")
plot1 <- VlnPlot(lung, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(lung, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
dev.off()

# Normalization, dimentional reduction, clustering and UMAP
lung <- SCTransform(lung, assay = "Spatial", verbose = FALSE)
lung <- RunPCA(lung, assay = "SCT", verbose = FALSE)
lung <- FindNeighbors(lung, reduction = "pca", dims = 1:30)
lung <- FindClusters(lung, verbose = FALSE)
lung <- RunUMAP(lung, reduction = "pca", dims = 1:30)

# UMAP & spatial plots
pdf("UMAP_cluster.pdf", width = 18, height = 8)
p1 <- DimPlot(lung, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(lung, label = TRUE, label.size = 3, repel = TRUE)
p1 + p2
dev.off()

# DEG analysis
lung.markers <- FindAllMarkers(lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- lung.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, "top10_visium.csv",append = FALSE, quote = TRUE, sep = "," , row.names = TRUE, col.names = NA)

# Save the objects
saveRDS(lung, file = "lung_visium.rds")

rm(list = ls())
gc();gc()
