#library(dplyr)
#library(reshape2)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(patchwork)
library(SPATA2)

set.seed(1234)
## Read data

Fun_lotate <- function( loc , do){
    c = do *pi/180
    x = loc[1]
    y = loc[2]
    x2 = x*cos(c) - y * sin(c)
    y2 = x*sin(c) + y * cos(c)
    return(c(x2,y2))
}

source("./Function-PlotTrraject.R")

obj <- "../demo_data/lung_visium.rds"
project = "LUAD_No3B"  # commandArgs(trailingOnly=TRUE)[1] # name
PAGE_res <- "../STEP1_Giotto_PAGE//test_PAGEscore.txt"
SeuratOBJ <- readRDS(obj)

# Read PAGE Result
PAGE_inf <- read.delim( PAGE_res  ,header = T,stringsAsFactors=F ,sep="\t")
rownames(PAGE_inf) <- PAGE_inf$cell_ID
PAGE_inf <- PAGE_inf[ , c(-1,-2) ]
# PAGE_inf$Annotation <- apply( PAGE_inf,1,function(x){  names(x)[which.max(x)] })

# Read Color Info
typeCol <- read.table("../demo_data//CelltypeCol.txt" , sep="\t" , header= T )
col_names <-  typeCol$col
names(col_names) <-  typeCol$CellType
Annotation<- factor(typeCol$CellType ,levels=typeCol$CellType)


bc <- intersect( rownames(PAGE_inf) , rownames(SeuratOBJ@meta.data))
bc_order <- rownames(SeuratOBJ@meta.data)[ rownames(SeuratOBJ@meta.data) %in% bc  ]
## Add PAGE score to Seurat OBJ
#for (i in colnames(PAGE_inf)){
for (i in c(1: (length(Annotation)  ) )  ) {
	Class =  as.character(Annotation[i])
    SeuratOBJ@meta.data[,Class] <- 0
    SeuratOBJ@meta.data[ bc_order ,Class ] <-  PAGE_inf[bc_order, Class]
}
SeuratOBJ@meta.data$Annotation <- PAGE_inf$CellType

# Convert SPATA OBJ
spata_obj <- transformSeuratToSpata( SeuratOBJ, project  , assay_name = 'Spatial', coords_from = "umap")
spata_obj <- adjustDirectoryInstructions( object = spata_obj,to = "spata_object", directory_new = "spata-obj-test.RDS" )

#Optional :  角度を調整する
coord_inf <- spata_obj@coordinates[[project]][,c(3,4)]
new_coord_inf <- t(apply( coord_inf  , 1 , Fun_lotate,do = -180  ))
new_coord_inf <-  t(apply(coord_inf ,1,  function(x){ return(x* c(1,-1))})) # y軸対象変換
spata_obj@coordinates[[project]][,c(3,4)] <- new_coord_inf

## make Cluster Centroid
center_list <- c()
all_clus <- unique( spata_obj@fdata[[project]][,"seurat_clusters"])
for ( i in all_clus  ) {
	cluster_cell <- spata_obj@coordinates[[project]][ spata_obj@fdata[[project]][,"seurat_clusters"] == i ,  ]
	center <-apply( cluster_cell[,c(3,4)],2,mean )
	center_id <-names(which.min(as.matrix(dist( rbind( center ,cluster_cell[,c(3,4)] ) ))[1,-1]))
	center_list <- c( center_list, center_id )
}
names(center_list) <- all_clus

#---------------------------------------------------------------------------------------------------
##
Tname <-"P1"
spata_obj <- createTrajectoryManually( spata_obj , trajectory_name = Tname ,
                         start = as.vector( new_coord_inf[ as.numeric(center_list["9"]), ]) ,
                         end   = as.vector( new_coord_inf[ as.numeric(center_list["3"]), ]) ,
                         vertices = list(
                                         v1 =  as.vector( new_coord_inf[ as.numeric(center_list["6"]), ])
                                         ),
                         width = 25,plot=T
                         )

p<-plotTrajectory(object = spata_obj,
               trajectory_name = Tname,
               color_by = "Annotation",
               pt_alpha = 0.25, # reduce alpha to highlight the trajectory's course
               display_image = FALSE) + legendTop() + ggplot2::scale_color_manual( values = col_names )
filePath <- paste0("./Spatial_Trajectory", project,"_", Tname ,".pdf")
ggsave(file = filePath, plot = p, dpi=100, width=8, height=8)

slope_min <- 0.04
Make_SCCR(spata_obj= spata_obj ,PAGE_inf=PAGE_inf ,project=project  ,trajectory_name=Tname ,
		 slope_min = slope_min , span =0.2 ,overlap = 4 , region_length =4 , Annotation = Annotation)

#--------------------------------------------------------
source("../GetTrajectory_Region.r")
ExtractBC <- getCellID( project,Tname,spata_obj  )
region_list <- unique(ExtractBC$region)
GeneExp <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
for(regionID in region_list ) {
    Spot_list <- ExtractBC %>% filter( region == regionID) %>% select(barcodes) %>% unlist
    Region_name <- ExtractBC %>% filter( region == regionID) %>% select(anno ,region) %>% distinct %>% unlist
    Region_name["region"] <- paste0( Tname ,"-", Region_name["region"] )

    tmpGeneExp  <- SeuratOBJ@assays$SCT@data[c("SPP1","FABP4") , Spot_list ] %>% apply(. , 1, mean)
    GeneExp <- rbind(GeneExp ,  c( project,"Invasive"  , Region_name , tmpGeneExp ) )

}

colnames( GeneExp) <-c("sample","Type" ,"region","regionID" ,"SPP1","FABP4")
write.table( GeneExp  , file=paste0("SPP1_FABP4_", project ,"_candidate_mat.txt"), quote=F, col.names=NA , sep="\t")

#-----------------------------------------------------------------------------
