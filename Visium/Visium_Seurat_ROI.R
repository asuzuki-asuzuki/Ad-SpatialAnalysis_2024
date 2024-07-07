library(Seurat) # v4.3.0
library(ggplot2)
library(patchwork)
library(dplyr)

library(RColorBrewer)
library(scales)

set.seed(1234)

lung <- readRDS("lung_visium.rds") # loading Seurat obj (LUAD No. 2 FFPE section C)

#####################################################################################################################
# Visualization of the boundary region between well-differentiated (NKX2-1+) and de-differentiated (HNF4A+) regions #
#####################################################################################################################

# Spatial plots with indicating the ROI in the dashed frame
pdf("zoom_1.pdf", width=14, height=7, pointsize = 18)
SpatialFeaturePlot(lung, features = c("NKX2-1", "HNF4A"), stroke = NA, image.alpha = 0, pt.size.factor = 2, ncol=2) & scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(0, NA), oob=squish) & NoGrid() & geom_rect(aes(xmin = 200, xmax = 400, ymin = 300, ymax = 450), fill = NA, color = "black", linetype="dashed", linewidth=1.5)
dev.off()

# Spatial plots (zoom-in the ROI)
max.y <- max(lung@images$slice1@coordinates$imagerow)*lung@images$slice1@scale.factors$lowres
min.y <- min(lung@images$slice1@coordinates$imagerow)*lung@images$slice1@scale.factors$lowres
ROI <- subset(lung, slice1_imagerow > max.y - 450 + min.y & slice1_imagerow < max.y - 300 + min.y & slice1_imagecol < 400 & slice1_imagecol > 200)
ratio = (450-300)/(400-200)

pdf("zoomin_1.pdf", width=21, height=14, pointsize = 18)
SpatialFeaturePlot(ROI, features = c("NKX2-1", "HNF4A", "IFITM1", "IFI6", "RHOV", "IDO1"), stroke = NA, image.alpha = 0, pt.size.factor = 3.5, ncol=3) & scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(0, NA), oob=squish) & NoGrid() & theme(aspect.ratio = ratio)
dev.off()

#######################################################################
# Visualization of the region from mucinous/proliferative to invasive #
#######################################################################

# Spatial plots with indicating the ROI in the dashed frame
pdf("zoom_2.pdf", width=14, height=7, pointsize = 18)
SpatialFeaturePlot(lung, features = c("MUC5AC", "SPINK1"), stroke = NA, image.alpha = 0, pt.size.factor = 2, ncol=2) & scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(0, NA), oob=squish) & NoGrid() & geom_rect(aes(xmin = 300, xmax = 480, ymin = 100, ymax = 320), fill = NA, color = "black", linetype="dashed", linewidth=1.5)
dev.off()

max.y <- max(lung@images$slice1@coordinates$imagerow)*lung@images$slice1@scale.factors$lowres
min.y <- min(lung@images$slice1@coordinates$imagerow)*lung@images$slice1@scale.factors$lowres
ROI <- subset(lung, slice1_imagerow > max.y - 320 + min.y & slice1_imagerow < max.y - 100 + min.y & slice1_imagecol < 480 & slice1_imagecol > 300)
ratio = (320-100)/(480-300)

# Spatial plots (zoom-in the ROI)
pdf("zoomin_2.pdf", width=21, height=14, pointsize = 18)
SpatialFeaturePlot(ROI, features = c("ACTA2", "MUC5AC", "MMP7", "CXCL14", "SPINK1"), stroke = NA, image.alpha = 0, pt.size.factor = 3.5, ncol=3) & scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(0, NA), oob=squish) & NoGrid() & theme(aspect.ratio = ratio, legend.position="right")
dev.off()


rm(list = ls())
gc();gc()
