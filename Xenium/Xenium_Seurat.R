library(Seurat) # v4.3.0
library(patchwork)
library(dplyr)

set.seed(1234)

# Please input the output directory of Xenium Analyzer
DATA <- "<INPUT>"

# Load output of Xenium Analayzer
#lung <- LoadXenium(DATA, fov = "fov")

assay <- 'Xenium'
fov <- 'fov'
data <- ReadXenium(
    data.dir = DATA,
    type = c("centroids", "segmentations"),
)
segmentations.data <- list(
    "centroids" = CreateCentroids(data$centroids),
    "segmentation" = CreateSegmentation(data$segmentations)
)
coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = assay
)
xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
xenium.obj[[fov]] <- coords
lung <- xenium.obj
rm(data)
rm(coords)
rm(segmentations.data)
rm(fov)
rm(xenium.obj)

lung <- subset(lung, subset = nCount_Xenium > 0)

# QC plot
pdf("Xenium_nCount.pdf")
VlnPlot(lung, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
dev.off()

# Normalization, dimentional reduction, clustering and UMAP
lung <- SCTransform(lung, assay = "Xenium")
lung <- RunPCA(lung, npcs = 30, features = rownames(lung))
lung <- RunUMAP(lung, dims = 1:30)
lung <- FindNeighbors(lung, reduction = "pca", dims = 1:30)
lung <- FindClusters(lung, resolution = 0.3)

# UMAP and spatial plots
pdf("UMAP_cluster.pdf", width = 18, height = 8)
p1 <- DimPlot(lung, reduction = "umap", label = TRUE)
p2 <- ImageDimPlot(lung, size = 0.3, border.size = NA, axes = TRUE) + NoGrid()
#p2 <- p2 & coord_flip() & theme(aspect.ratio = max(p2$data$y)/max(p2$data$x)) & scale_x_reverse() & scale_y_reverse() & ggtitle(NULL) # flipping if needed
p1 + p2
dev.off()

# Save the object
saveRDS(lung, file = "lung_xenium.rds")

## DEG analysis
lung.markers <- FindAllMarkers(lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- lung.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, "top10_xenium.csv",append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)

# Save the object
saveRDS(lung.markers, file = "lung_marker.rds")

rm(list = ls())
gc();gc()
