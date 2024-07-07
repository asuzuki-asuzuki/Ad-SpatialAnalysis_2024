library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(reshape2)
library(cowplot)
library(patchwork)

library(tidyverse)
library(skimr)
library(summarytools)
library(psych)
library(pheatmap)

## Read data
WHO_GEX.list <- NULL

for (i in list.files("./data/" ) ) {
        WHO_GEX.list <- c(WHO_GEX.list, readRDS( paste0("./data/",i) ) )
}

names(WHO_GEX.list) <-  c("FFPE_2A","FFPE_2B","FFPE_2C","FFPE_2D")
FFPE.Merge <- merge(x = WHO_GEX.list[[1]], y = unlist(WHO_GEX.list)[2:4],add.cell.ids =names(WHO_GEX.list) , project = "Visium")

DefaultAssay(FFPE.Merge) <- "Spatial"
FFPE.Merge[["sample"]] <- paste( c(FFPE.Merge@meta.data$Protcol) ,  c(FFPE.Merge@meta.data$Area) , sep="_" )

project<-"FFPE_merge"

filePath <- paste0("./Raw_Countplot_", project ,".pdf")
p<-VlnPlot( FFPE.Merge , features = c("nCount_Spatial"), group.by="sample"  , ncol = 1, pt.size = 0) + NoLegend()
ggsave(file = filePath, plot = p, dpi=100, width=8, height=8)
print(dim(FFPE.Merge))

filePath <- paste0("./Raw_Feature_plot_", project ,".pdf")
p1<- VlnPlot( FFPE.Merge , features = c("nFeature_Spatial"), ncol = 1, pt.size = 0) + NoLegend()
p2<- SpatialFeaturePlot (  FFPE.Merge , features = c("nFeature_Spatial" )) + theme(legend.position = "right")
p <- wrap_plots(p1, p2)
ggsave(file = filePath, plot = p, dpi=100, width=8, height=8)

FFPE.Merge@meta.data %>% group_by(sample) %>% summarise( mean_geture = median(nFeature_Spatial))


DefaultAssay(FFPE.Merge) <- "Spatial"
FFPE.Merge <- FindVariableFeatures(FFPE.Merge, selection.method = "vst", nfeatures =4000)
all.genes <- rownames(FFPE.Merge)
FFPE.Merge <- ScaleData(FFPE.Merge , features = all.genes , assay = "Spatial", verbose = FALSE)
FFPE.Merge <- RunPCA(FFPE.Merge, assay = "Spatial", verbose = FALSE)

filePath="./VizDimLoadings.png"
p<-VizDimLoadings(FFPE.Merge, dims = 1:5, reduction = "pca")
ggsave(file = filePath, plot = p, dpi=50, width=16, height=16)
filePath="./PCAPlot.png"
p<-DimPlot(FFPE.Merge, reduction = "pca", group.by="orig.ident")
ggsave(file = filePath, plot = p, dpi=100, width=16, height=16)
filePath="./PCheatmap.png"
p<-DimHeatmap(FFPE.Merge, dims = 1:15, cells = 500, balanced = TRUE,fast = FALSE)
ggsave(file = filePath, plot = p, dpi=100, width=24, height=24)
filePath="./ElbowPlot.png"
p<-ElbowPlot(FFPE.Merge)
ggsave(file = filePath, plot = p, dpi=100, width=16, height=8)


pcadim=c(1:30)
resol = 0.5
FFPE.Merge <- FindNeighbors(FFPE.Merge, reduction = "pca", dims = pcadim)
FFPE.Merge <- FindClusters(FFPE.Merge, resolution = resol)
FFPE.Merge <- RunUMAP(FFPE.Merge, reduction = "pca", dims = pcadim)

freq_table <- prop.table(x = table(Idents(FFPE.Merge), FFPE.Merge@meta.data[, "sample"]), margin = 2)
freq_mmod<-melt(t(freq_table),id=c("sample","cluster","val"))
colnames(freq_mmod) <- c("sample","cluster","val")
freq_mmod$cluster <- as.factor( freq_mmod$cluster )

sink("frek_table_RNA.txt")
table(Idents(FFPE.Merge), FFPE.Merge@meta.data[, "orig.ident"])
sink()

sink("freq_tableprop_RNA.txt")
freq_table
sink()

filePath="./umap_RNA_seuraCluster.pdf"
p1 <- DimPlot(object = FFPE.Merge , group.by="seurat_clusters",reduction = "umap", pt.size = 0.2 ,
        label=T )+ ggtitle(label = "UMAP cluster")
p2 <- DimPlot(object = FFPE.Merge , group.by="sample",reduction = "umap", pt.size = 0.2 ,
		label=T )+ ggtitle(label = "Sample")
p <- p1 + p2
ggsave(file = filePath , plot = p, dpi=100, width=12, height=6 )

Idents(FFPE.Merge) <- "seurat_clusters"
filePath="./Spatial_plot.pdf"
p2 <- SpatialDimPlot( FFPE.Merge ,label = TRUE, label.size = 10,ncol=2)
ggsave(file = filePath , plot = p2, dpi=100, width=25, height=25 )


image_list <- names(FFPE.Merge@images)
name_list  <-  c("FFPE_2A","FFPE_2B","FFPE_2C","FFPE_2D")

for( i in c(1:length(image_list))){
 p<- SpatialDimPlot( FFPE.Merge ,label = F, images =image_list[i] , pt.size.factor =2) # + ggtitle(name_list[i])
 ggsave(file = paste0("./Spatial_",name_list[i],"_plot.pdf"),plot = p, dpi=100, width=5, height=5 )

}

# Scatterig Plot
avg.sample <- log(as.data.frame(AverageExpression(FFPE.Merge, verbose = FALSE , assays = "Spatial",group.by = "sample"))+1,2)

pdf("./FFPE_scatterplot.pdf", width=8, height=8)
pairs.panels( avg.sample ,hist.col="white",rug=F,ellipses=F,lm=T,cex.cor =1,cex=0.5)
dev.off()


Idents(FFPE.Merge) <- "seurat_clusters"
for (i in  unique(Idents(FFPE.Merge))){
	print(i)
	cluster.cells <- subset(FFPE.Merge, idents = i)
	cluster.avg.sample <- log(as.data.frame(AverageExpression(cluster.cells, verbose = FALSE ,
		assays = "Spatial",group.by = "sample"))+1,2)

	pdf( paste0("./FFPE_Clustr",i,"_scatterplot.pdf") , width=8, height=8)
	print(pairs.panels( cluster.avg.sample ,hist.col="white",rug=F,ellipses=F,lm=T,cex.cor =1,cex=0.5))
	dev.off()	
}


saveRDS( FFPE.Merge ,file="FFPE_Data.rds")



MerkerGene <- read.table("./marker_mini.txt" , sep="\t" , header= T )
MyPlot <- function( OBJ , GeneSet , prj = "prg" , pmethod ="umap",cluster = "seurat_clusters" , width=6, height=6 ){
    filePath = paste0("./",prj,"_gene_Plot_EXPplot",".pdf")
    col_num = length(GeneSet)

    p<-VlnPlot(object = OBJ, features=GeneSet, pt.size = 0, ncol=col_num, group.by=cluster, split.by="sample")
    ggsave(file = filePath, plot = p, dpi=100, width=width, height=height)
    print(paste0("finish ", prj ," plot!!" ) )
}

for( i in unique(MerkerGene$Celltype)){
	g_list <-  MerkerGene[MerkerGene[,2] == i, 1]
	MyPlot( OBJ = FFPE.Merge, GeneSet= g_list, prj =i, width=20, height=4)
}

p<-VlnPlot(object = FFPE.Merge, features=MerkerGene[,1], pt.size = 0, ncol=4, split.by="sample")
ggsave(file = "./MarkerGeneEXPplot.pdf", plot = p, dpi=100, width=10, height=40)


SpatialFeaturePlot(FFPE, features = c("CD68","SFTPB","MUC5AC","SPINK1", alpha = c(0.1, 1)))
SpatialFeaturePlot(FFPE, features = c("IGKC", alpha = c(0.1, 1)))

p<-VlnPlot(object =  FFPE.Merge , features= c("CD68","SFTPB","MUC5AC","SPINK1","SOD2"), pt.size = 0,ncol=2 ,
         group.by="seurat_clusters",split.by="sample" )
filePath="./FFPE_MarkerVinplot.pdf"
ggsave(file = filePath , plot = p, dpi=100, width=10, height=10 )


image_list <- names(FFPE.Merge@images)
name_list  <-  c("FFPE_2A","FFPE_2B","FFPE_2C","FFPE_2D")

for( i in c(1:length(image_list))){
	for( g in c("CD68","SFTPB","MUC5AC","SPINK1","SOD2") ) {
	 	p<- SpatialFeaturePlot( FFPE.Merge , features = g , images =image_list[i] ) # + ggtitle(name_list[i])
		ggsave(file = paste0("./SpatialfeturePlot_",g,"_",name_list[i],"_plot.pdf"),plot = p, dpi=100, width=5, height=5 )
	}
}


    for( g in c("CD68","SFTPB","MUC5AC","SPINK1","SOD2") ) {
        p<- SpatialFeaturePlot( FFPE.Merge , features = g  ) # + ggtitle(name_list[i])
        ggsave(file = paste0("./SpatialfeturePlot_4sample_",g,"_plot.pdf"),plot = p, dpi=100, width=20, height=5 )
    }


#---- 

library(pheatmap)
FFPE.Merge <-readRDS("FFPE_Data.rds")

FFPE.Merge[["Sample_Cluster"]] <- paste( FFPE.Merge@meta.data$sample, FFPE.Merge@meta.data$seurat_clusters )
ave_exp <-log(as.data.frame(AverageExpression(FFPE.Merge, verbose =FALSE, assays="Spatial",group.by="Sample_Cluster")) +1 ,2 )

write.table(round(cor(ave_exp),2) , file="SampleClusterCor.txt" , sep="\t" , col.names=NA , row.names=T , quote=F)

cor_data <-  cor(ave_exp)
cluste_inf <- as.data.frame( cbind(sample = gsub("Spatial.|.[0-9]+","", rownames(cor_data)) , 
				cluster=gsub("Spatial.FFPE_[ABCD].","", rownames(cor_data)) ))
rownames( cluste_inf ) <-  rownames(cor_data)


new_col <- gsub("Spatial.","",colnames(cor_data))
new_col_order<-paste(  rep(c("FFPE_A","FFPE_B","FFPE_C","FFPE_D"),13) ,rep(c(0:13),each=4) ,sep="." )

colnames(cor_data) <- new_col
rownames(cor_data) <- new_col

cor_data <- cor_data[new_col_order[-49],new_col_order[-49]]


pdf( paste0("./FFPE_Sample_Clustr_corheatmap.pdf") , width=10, height=10)
pheatmap( cor_data , scale = "none" ,     
        show_rownames =T , show_colnames= T, cluster_rows = F,cluster_cols = F
        ) 
dev.off()


cor_data2 <- cor_data
for (i in c(1:54 ) ){
	cor_data2[i,(i+1):55] <- NA
}

pdf( paste0("./FFPE_Sample_Clustr_corheatmap2.pdf") , width=10, height=10)
pheatmap( cor_data2 , scale = "none" , 
        show_rownames =T , show_colnames= T, cluster_rows = F,cluster_cols = F
        ) 
dev.off()

