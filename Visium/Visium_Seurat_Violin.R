library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

lung <- readRDS("lung_visium.rds") # Visium Seurat obj

pdf("violin.pdf", width = 18, height = 16)
VlnPlot(lung, features = c("NAPSA", "SFTPB", "DUOX1", "IFITM1", "MUC5AC", "SLC2A1", "SPINK1", "MMP7", "CXCL14", "ACTA2", "SPARC", "COL1A1", "COL1A2"), ncol=3)
dev.off()

rm(list = ls())
gc();gc()
