library(dplyr)
library(tidyverse)
library(ggplot2)
library(gggibbous)
library(sf)

#getCellID( project,trajectory_name,spata_obj  )
extrac_cate = c( "Proliferative","Invasive" )
getCellID<- function( project,trajectory_name,spata_obj , extrac_cate ){
	#print(trajectory_name)
	filename = paste0( "ScoreChangedRegoin_", project,"_", trajectory_name ,".txt")
	SCCR_list <- read.delim( file= filename, header = T,stringsAsFactors=F ,sep="\t",row.name=1)  %>%
		filter( Annotaion %in% extrac_cate  ) %>%
		select(x_start , x_end , pos , Annotaion  ) %>%
		distinct()

	tjkMat    <- read.delim( file= paste0( "trajectory_score_matrix_", project,"_", trajectory_name ,".txt") , header = T,stringsAsFactors=F ,sep="\t",row.name=1)  
	segment_inf <-  spata_obj@trajectories[[project]][[trajectory_name]]@segment_trajectory_df
	segment_inf$part<- gsub( "_" , " " , gsub("p", "P", segment_inf$part))
	segment_inf$length <-  apply(segment_inf, 1, function(x , tjkMat){ return( sum(tjkMat$trajectory_part == x["part"]))},tjkMat  )
	segment_inf <- cbind(segment_inf, t(apply(segment_inf,1,function(x){ 
		x=as.double( x[c("x","y","xend","yend","length")]);return(c( xbin =(x[3]-x[1])/x[5] , ybin=(x[4]-x[2])/x[5] ))
	}) ) )

	segment_inf$st <- tjkMat%>% select(trajectory_part,trajectory_part_order, trajectory_order) %>% 
		filter( trajectory_part_order == 1 ) %>% select(trajectory_order  ) %>% unlist
	segment_inf$en  <- segment_inf$st + segment_inf$length -1

	# checke belongs Part :
	aria_bin <- 40
	region_id <- 0
	ExtractBarcode <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]

	for ( SCCR_id in c(1:dim(SCCR_list)[1]  ) ) {
		region_id <- region_id +1
		SCCR <-  SCCR_list[SCCR_id,]
		seurat_meta <- spata_obj@trajectories[[project]][[trajectory_name]]@compiled_trajectory_df
		seurat_meta$region <- "NA"
		seurat_meta$anno   <- SCCR$Annotaion

		for( part_id in c(1:dim(segment_inf)[1])  ){
			part <- segment_inf[part_id ,]
			new_part <- part
			if(  SCCR$x_start  < part$en & part$st < SCCR$x_end  ) {
				# new_part
				if( part$st < SCCR$x_start){
					new_part$x <- part$x  + ( part$xbin * ( SCCR$x_start - part$st ))
					new_part$y <- part$y  + ( part$ybin * ( SCCR$x_start - part$st ))
					new_part$st <-  SCCR$x_start
				} 
				if ( SCCR$x_end  < part$en ) {
					new_part$xend <- part$x  + ( part$xbin * ( part$en - SCCR$x_start ))
					new_part$yend <- part$y  + ( part$ybin * ( part$en - SCCR$x_end ))
					new_part$en   <-  SCCR$x_end
				}
				print( new_part  )
				# 直行専
				slope     <- (new_part$yend -  new_part$y) / ( new_part$xend - new_part$x )
				intercept <- new_part$y - (slope * new_part$x)
				line_eq <- sprintf("y = %.4fx + %.4f", slope, intercept)

				perpendicular_slope <- -1 / slope
				perpendicular_intercept_st <-  new_part$y    - (perpendicular_slope * new_part$x)
				perpendicular_intercept_en <-  new_part$yend - (perpendicular_slope * new_part$xend)

				zone <- list (  A= c(x = new_part$x - aria_bin, y = perpendicular_slope* ( new_part$x - aria_bin ) + perpendicular_intercept_st ) ,  # 左下
								B= c(x = new_part$x + aria_bin, y = perpendicular_slope* ( new_part$x + aria_bin ) + perpendicular_intercept_st ) ,  # 左上
								C= c(x = new_part$xend - aria_bin, y = perpendicular_slope* ( new_part$xend - aria_bin ) + perpendicular_intercept_en ),  # 右下
								D= c(x = new_part$xend + aria_bin, y = perpendicular_slope* ( new_part$xend + aria_bin ) + perpendicular_intercept_en ) )  # 右上

				polygon_dt <- st_polygon( list(rbind( zone[["A"]] , zone[["B"]] , zone[["D"]],zone[["C"]] ,zone[["A"]] ) ) )
				polygon_sf <- st_sf(geometry = st_sfc(polygon_dt))
				point_df <- seurat_meta %>% filter( trajectory_part == new_part$part   ) %>% select( barcodes , x ,y )
				
				points_sf <- st_as_sf(point_df, coords = c("x", "y"), crs = st_crs(polygon_sf))
				inside_points <- st_intersection(points_sf, polygon_sf) %>% as.data.frame %>% select(barcodes)

				seurat_meta [seurat_meta$barcodes %in% inside_points$barcodes ,"region"] <- paste0( "C", region_id)
				
			}
		}
		seurat_meta %>% filter( region != "NA") %>% select(barcodes , region ,anno)
		ExtractBarcode <-rbind(ExtractBarcode ,  seurat_meta %>% filter( region != "NA") %>% select(barcodes , region ,anno) %>% as.data.frame )
	}
	return(  ExtractBarcode)
}


