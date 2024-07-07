library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dbscan)
library(tidyr)
library(tidyverse)
library(ggridges)
library(sf)
set.seed(1234)

source("Function_plotModule.r")

#++++++++++++++++++++++++++++++++++++++++++++++
## Read Input File
#++++++++++++++++++++++++++++++++++++++++++++++

Sample <- "test" 
obj           <- paste0("../demo_data/lung_visium.rds" )
PAGE_inf.file <- paste0("../STEP1_Giotto_PAGE/test_PAGEscore.txt")
typeCol       <- read.table("../demo_data/CelltypeCol.txt" , sep="\t" , header= T )
SeuratOBJ     <- readRDS(obj) 
PAGE_inf      <- read.table(PAGE_inf.file , sep="\t" , header= T , row.names=1)
rownames(PAGE_inf) <- PAGE_inf$cell_ID

# spot coordinates matrix
mat.dt     <- SeuratOBJ@images$slice1@coordinates %>%
				select( c(imagecol,imagerow)) %>% 
				as.data.frame
mat.dt     <- mat.dt[ rownames(PAGE_inf) , ]
mat.dt[,2] <- mat.dt[,2]*(-1) # conver angel
colnames(mat.dt) <- c("x","y")
mat.dt <-cbind(mat.dt, PAGE_inf) # add TME score
PAGE_inf <- PAGE_inf[ , c(2:(ncol(PAGE_inf)-1 ) ) ]

col_names        <- typeCol$col
names(col_names) <- typeCol$CellType

#++++++++++++++++++++++++++++++++++++++++++++++
# Density Score
#++++++++++++++++++++++++++++++++++++++++++++++
# spot毎の距離重み付け集計計算#################################################################
# spot間のpixel距離を求める

forCheckSpotDist <- frNN(x= mat.dt[,c("x","y")], eps =2000)
summary(unlist( lapply(forCheckSpotDist$dist,min) ))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  131.0   131.0   131.0   131.0   131.0   228.1

max_step =35
tmp_eps <- 134 # 要確認
cate_num <- length( names(col_names))

# This process takes a few hours.
nnCounts_list  <-c()
for (size_c in c(1:max_step)){
    eps <- tmp_eps * size_c
    nn_data <- frNN(x= mat.dt[,c("x","y")], eps = eps)
    tmp_data<- lapply( 
				rownames(mat.dt) , 
				RegionScore , 
				nn_data = nn_data ,
				mat.dt =mat.dt,
				tmp_eps = tmp_eps,
				size_c= size_c ,
				col_names = col_names  
				)
    nnCounts_list  <- rbind(nnCounts_list , do.call(rbind ,tmp_data ) )
}

nnCounts_list <- nnCounts_list %>%
                 mutate( 
					x = mat.dt[ID , "x" ] , 
					y = mat.dt[ID , "y" ] , 
					Annotation = mat.dt[ID , "CellType" ] 
				 )

########################################################
## Define Enrichment Spot
########################################################
weight= 0.7
lowlim = 10

dot_size=1.2
RegionMinSize = 15


PlotAnnotation <- "Proliferative"
max_score <- nnCounts_list %>% 
			 filter(CellType == PlotAnnotation ) %>%
			 select(adj_counts) %>% 
			 max

threshold = ifelse( 
				max_score * weight > lowlim ,
				round(max_score * weight, 2) ,
				lowlim 
			)

Aggressive_plot_res <- makeRegionPlot( 
							nnCounts_list= nnCounts_list ,
							col_names=col_names,
							Plot_annotation=PlotAnnotation,
							dot_size= 2.5,
							p1_lowCol= "gray",
							p1_highCol="green",
							threshold = threshold,
							max_step=35 ,
							name = Sample ,
							width=22 ,
							height=7 ,
							dpi=200)

# Check the result of the previous makeRegionPlot() and if multiple regions are detected,
# specify each region manually by enclosing it with a rectangle (specifying the coordinates of the four corners)
# in order to cut out the relevant region.

# A = left_bottom ,B=left_top , C= right_bottom, D=rite_top
# A != D , B != C
frame_list <- list(
		        list ( A= c(x= 1100, y= -8500), B= c(x = 1100 , y= -6200), C=c(x= 4900 , y= -8500), D=c (x = 4900 , y=-6200) )
			)


Split_Res1<-PlotSplitRegion(
                frame_list=frame_list ,
                plot_res= Aggressive_plot_res ,
                name=Sample ,
                dot_size= dot_size,
                RegionMinSize = RegionMinSize ,
                Plot_annotation = "Aggressive",
                width=17,
                height=5,
                dpi=150 ,
                p1_lowCol="gray",
                p1_highCol="green" ,
                col_names = col_names
            )




PlotAnnotation <- "Invasive"
max_score <- nnCounts_list %>% 
			 filter(CellType == PlotAnnotation ) %>% 
			 select(adj_counts) %>%
			 max

threshold = ifelse( 
				max_score * weight > lowlim ,
				round(max_score * weight, 2) ,
				lowlim  
			)

Invation_plot_res <- makeRegionPlot( 
							nnCounts_list= nnCounts_list ,
							col_names=col_names,
							Plot_annotation=PlotAnnotation,
							dot_size= 2.5,
							p1_lowCol= "gray",
							p1_highCol="green",
							threshold = threshold,
							max_step=35 ,
							name = Sample ,
							width=22 ,
							height=7 ,
							dpi=200)


frame_list <-list(
        list ( A= c(x= 4500, y= -5400), B= c(x = 4500 , y= -4500), C=c(x= 5600 , y= -5400), D=c (x = 5600 , y=-4500) ),
        list ( A= c(x= 8000, y= -6400), B= c(x = 8000 , y= -3800), C=c(x= 9800 , y= -6400), D=c (x = 9800 , y=-3800) )
)

Split_Res2<-PlotSplitRegion(
                frame_list=frame_list ,
                plot_res= Invation_plot_res ,
                name=Sample ,
                dot_size= dot_size,
                RegionMinSize = RegionMinSize ,
                Plot_annotation = PlotAnnotation,
                width=17,
                height=5,
                dpi=150 ,
                p1_lowCol="gray",
                p1_highCol="green" ,
                col_names = col_names
            )







