library(Seurat) # v5.1.0
library(patchwork)
library(dplyr)

set.seed(1234)

# Please input the segmentation result of PhenoCyler (CSV format)
codex.obj <- LoadAkoya(filename = "<INPUT>", type = "qupath", fov = "phenocycler")

# number of cells
length(codex.obj$nCount_Akoya)

# Filtering: Sum of intensity: <99 percentile & >1 percentile
high <- quantile(codex.obj$nCount_Akoya, prob=0.99)
low <- quantile(codex.obj$nCount_Akoya, prob=0.01)
codex.obj <- subset(codex.obj, subset = nCount_Akoya > low  & nCount_Akoya < high)
# number of cells
length(codex.obj$nCount_Akoya)

# Normalization and scaling
codex.obj <- NormalizeData(object = codex.obj, normalization.method = "CLR", margin = 2)
codex.obj <- ScaleData(codex.obj)

# Setting variable features
VariableFeatures(codex.obj) <- rownames(codex.obj)

# Dimentional reduction, clustering and UMAP
codex.obj <- RunPCA(object = codex.obj, npcs = 10, verbose = FALSE)
codex.obj <- RunUMAP(object = codex.obj, dims = 1:10, verbose = FALSE)
codex.obj <- FindNeighbors(object = codex.obj, dims = 1:10, verbose = FALSE)
codex.obj <- FindClusters(object = codex.obj, verbose = FALSE, resolution = 0.4)

# UMAP & spatial plots
pdf("UMAP_cluster.pdf", width = 18, height = 8)
p1 <- DimPlot(codex.obj, reduction = "umap", label = TRUE)
p2 <- ImageDimPlot(codex.obj, size = 0.3, border.size = NA, axes = TRUE) + NoGrid()
#p2 <- p2 & coord_flip() & theme(aspect.ratio = max(p2$data$y)/max(p2$data$x)) & scale_x_reverse() & scale_y_reverse() & ggtitle(NULL) # flipping if needed
p1 + p2
dev.off()

# Save the object
saveRDS(codex.obj, file = "lung_phenocycler.rds")

# DEG analysis
lung.markers <- FindAllMarkers(codex.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- lung.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, "top10_phenocycler.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)

# Save the object
saveRDS(lung.markers, file = "lung_marker.rds")

rm(list = ls())
gc();gc()
