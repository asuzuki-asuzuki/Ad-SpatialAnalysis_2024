library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tidyverse)
library(dbscan)

Fun_lotate <- function( loc , do){
    c = do *pi/180
    x = loc[1]
    y = loc[2]
    x2 = x*cos(c) - y * sin(c)
    y2 = x*sin(c) + y * cos(c)
    return(c(x2,y2))
}


RegionScore <- function(id_name , nn_data , mat.dt , tmp_eps , size_c , col_names ){
    countdt <- cbind(
					dist= nn_data$dist[[id_name]] ,
					pos = nn_data$id[[id_name]]   , 
					mat.dt[ nn_data$id[[id_name]] , ] 
				) %>%
		        rbind( . , 
					cbind( 
						dist = tmp_eps ,
						pos = which( rownames(mat.dt) == id_name),
						mat.dt[ id_name , ] 
					) 
				)  %>% 
				as.data.frame() %>%
		        mutate( num =1 ,adj_num= num * tmp_eps/dist ) %>%
				group_by(CellType) %>% 
				summarise('counts' = sum(num) , 'adj_counts'=sum(adj_num))

    na_data <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]

    for ( chechname in names(col_names) ) {
        if( ! chechname %in% countdt$CellType ) {
            na_data <- rbind( na_data, c( chechname, as.double(0) , as.double(0) ) )
        }
    }

    colnames(na_data)  <- colnames(countdt)
    na_data$counts     <- as.double(na_data$counts)
    na_data$adj_counts <- as.double(na_data$adj_counts)
    if ( nrow(na_data) > 0 ){
        countdt <-  bind_rows( countdt ,  na_data  )
    }
    countdt %>%
	    mutate( 
			ID = id_name , 
			total = sum(counts) ,
			adj_total = sum(adj_counts) ,
			rate = counts/total,
			adj_rate = adj_counts/adj_total,
			step = size_c 
		) %>%  return()
}


#############################
# Detect PM region
############################


makeRegionPlot <- function(
						nnCounts_list ,
						col_names,
						Plot_annotation ,
						dot_size=2 ,
						p1_lowCol="gray",
						p1_highCol="green",
            			threshold =15 ,
						max_step=35 ,
						name="Sample", 
						width=17,
						height=5,
						dpi=100
){
    max_score_ID = nnCounts_list %>% 
							filter( step == max_step  &  CellType == Plot_annotation ) %>%
							filter( adj_counts == max(adj_counts) ) %>%
							select(ID) %>% as.vector
	# Plot local Enrichment score
    p1 <- nnCounts_list %>% 
			filter( CellType == Plot_annotation & step == max_step ) %>% 
			ggplot( aes( x=x ,y=y , color = adj_counts )) +
            geom_point( size = dot_size ) + 
			scale_color_gradient(low = p1_lowCol, high = p1_highCol) +
            theme(legend.position = "right", axis.title.x  =element_blank() ,axis.title.y  =element_blank()) +
            labs(title = paste0(Plot_annotation) )

	## Define a Enrichment spot as a spot where the Enrichment score exceeds a threshold value.
	## The spot with the highest SCORE is also defined.
    nnCounts_Region_add2bit <- nnCounts_list %>%  
								filter( step == max_step ,  CellType == Plot_annotation) %>%
								mutate( 
										Region = ifelse( ( adj_counts > threshold )  , "Target" ,"Other"  ),
										CenterSpot =  ifelse( ( ID == max_score_ID  ) , "Center" ,"Other" )
								)

	# Plot Enrichment Spot
    p2 <- nnCounts_Region_add2bit %>%
				ggplot( aes( x=x ,y=y , color = Region )) +
				geom_point(size = dot_size) + 
				scale_color_manual( values = c( "Target"= p1_highCol , "Other" = p1_lowCol ) ) +
				labs(title =  paste0(Plot_annotation, " SCore FIlter > " , threshold) )

	# 
    p3 <- nnCounts_Region_add2bit %>% 
			ggplot( aes( x=x ,y=y ) )  + 
			geom_point( aes( color=Region  , fill=  Annotation ), shape=21, size = dot_size ) +
            scale_fill_manual( values = col_names) + 
			scale_color_manual( values = c( "Target"= "red" , "Other" =  "white") )

    p<- p1+p2+p3
    filePath <- paste0( "./Sub_DensityScorePlot_",Plot_annotation,"_",name,".pdf")
    ggsave(file = filePath, plot = p , dpi=dpi, width=width, height=height)
    return(nnCounts_Region_add2bit)

}



PlotSplitRegion <-function( 
							frame_list ,
							plot_res ,
							name ,
							dot_size=2,
							RegionMinSize = 10 ,
							Plot_annotation , 
							width=17,
							height=5,
							dpi=100 ,
							p1_lowCol="gray",
							p1_highCol="green" ,
							col_names 
){
	plot_res$PlotRegion <- plot_res$Region
    modRegtionPlot <-plot_res %>%
						ggplot( aes( x=x ,y=y , color = Region )) +
						geom_point(size = dot_size)
	
	# A = left_bottom ,B=left_top , C= right_bottom, D=rite_top
	# A != D , B != C
    for (i in c(1:length(frame_list)) ) {
        frame_pos <- frame_list[[i]]
        modRegtionPlot <- modRegtionPlot +
            geom_segment(  x= frame_pos[["A"]]["x"] ,xend= frame_pos[["C"]]["x"] ,y= frame_pos[["A"]]["y"] ,yend=  frame_pos[["C"]]["y"] , color = "yellow" ) +
            geom_segment(  x= frame_pos[["B"]]["x"] ,xend= frame_pos[["D"]]["x"] ,y= frame_pos[["B"]]["y"] ,yend=  frame_pos[["D"]]["y"] , color = "yellow" ) +
            geom_segment(  x= frame_pos[["A"]]["x"] ,xend= frame_pos[["B"]]["x"] ,y= frame_pos[["A"]]["y"] ,yend=  frame_pos[["B"]]["y"] , color = "yellow" ) +
            geom_segment(  x= frame_pos[["C"]]["x"] ,xend= frame_pos[["D"]]["x"] ,y= frame_pos[["C"]]["y"] ,yend=  frame_pos[["D"]]["y"] , color = "yellow" ) 
    }

    # Filter
    plot_resSplit <- plot_res
    plot_resSplit$CenterID <- "Other"

    skip_list <-c()
    printCount = 0
    for (i in c(1:length(frame_list)) ) {
        frame_pos <- frame_list[[i]]
		# extract dot
		polygon_dt <- st_polygon( list(rbind( frame_pos[["A"]] , frame_pos[["B"]] , frame_pos[["D"]],frame_pos[["C"]] , frame_pos[["A"]] ) ) )
		polygon_sf <- st_sf(geometry = st_sfc(polygon_dt))
		point_df <- plot_resSplit %>% select(c(ID,x,y) ) %>% as.data.frame
		points_sf <- st_as_sf(point_df, coords = c("x", "y"), crs = st_crs(polygon_sf))
		inside_points <- st_intersection(points_sf, polygon_sf) %>% as.data.frame %>% select(ID)
        Regionsize <- plot_resSplit %>% filter( ID %in% inside_points$ID & Region == "Target") %>% nrow

		plot_resSplit

        if (Regionsize <= RegionMinSize ) {
            plot_resSplit <- plot_resSplit %>%
				mutate(PlotRegion = if_else(  PlotRegion == "Target" & ID %in% inside_points$ID , "Other" , PlotRegion ) ) %>%
				mutate(Region = if_else(  Region == "Target" & ID %in% inside_points$ID , "Other" , Region ) )

                skip_list <- c(skip_list , i)
        }else {
            printCount =printCount + 1
			center_id <- plot_resSplit %>% filter(ID %in% inside_points$ID) %>% filter( adj_counts == max( adj_counts ))%>% select(ID) %>% unlist

            plot_resSplit <- plot_resSplit %>%
				mutate(PlotRegion = if_else(  PlotRegion == "Target" & ID %in% inside_points$ID , "TargetRegion" , PlotRegion ) ) %>% 
				mutate(CenterID   = if_else( ID ==center_id   , paste0("Center-",printCount) , CenterID  ) ) %>% 
				mutate(Region     = if_else(  Region == "Target" & ID %in% inside_points$ID , paste0("Region-",printCount) , Region ) )
        }
    }

	plot_resSplit <- plot_resSplit %>% mutate(Region = if_else(  Region == "Target" , "Other"  , Region ) , PlotRegion = if_else( PlotRegion == "Target"  , "Other"  , PlotRegion  ) )

    modRegtionPlot2 <- plot_resSplit %>%ggplot( aes( x=x ,y=y , color = PlotRegion )) + geom_point(size = dot_size)
    printCount = 0
    for (i in c(1:length(frame_list)) ) {
        if ( ! (i %in% skip_list) ){
            printCount =printCount + 1
            frame_pos <- frame_list[[i]]
            modRegtionPlot2 <- modRegtionPlot2 +
	            geom_segment(  x= frame_pos[["A"]]["x"] ,xend= frame_pos[["C"]]["x"] ,y= frame_pos[["A"]]["y"] ,yend=  frame_pos[["C"]]["y"] , color = "yellow" ) +
    	        geom_segment(  x= frame_pos[["B"]]["x"] ,xend= frame_pos[["D"]]["x"] ,y= frame_pos[["B"]]["y"] ,yend=  frame_pos[["D"]]["y"] , color = "yellow" ) +
        	    geom_segment(  x= frame_pos[["A"]]["x"] ,xend= frame_pos[["B"]]["x"] ,y= frame_pos[["A"]]["y"] ,yend=  frame_pos[["B"]]["y"] , color = "yellow" ) +
            	geom_segment(  x= frame_pos[["C"]]["x"] ,xend= frame_pos[["D"]]["x"] ,y= frame_pos[["C"]]["y"] ,yend=  frame_pos[["D"]]["y"] , color = "yellow" ) +
                annotate("text", x= frame_pos[["B"]]["x"]+200 ,   y= frame_pos[["B"]]["y"] +200 , label= paste0("Region-",printCount))
        }
    }

    p1 <- plot_res %>% ggplot( aes( x=x ,y=y , color = adj_counts )) +
            geom_point( size = dot_size ) + scale_color_gradient(low = p1_lowCol, high = p1_highCol) +
            theme(legend.position = "right", axis.title.x  =element_blank() ,axis.title.y  =element_blank()) +
            labs(title = "Neighbor score" , subtitle = paste0(Plot_annotation) )

	modRegtionPlot  <- modRegtionPlot  +  labs(title = "Befor Filter" , subtitle = Plot_annotation)
	modRegtionPlot2 <- modRegtionPlot2 +  labs(title = paste0("After Filter: Region Size <= " , RegionMinSize ) , subtitle = Plot_annotation)
	p<-p1 + modRegtionPlot + modRegtionPlot2
    filePath <- paste0( "./RegionFilterPlot_",Plot_annotation,"_",name,".pdf")
    ggsave(file = filePath, plot = p , dpi=dpi, width=width, height=height)


	PopulationInRegion <- prop.table( table(plot_resSplit$Region , plot_resSplit$Annotation  ) ) %>% as.data.frame %>% 
					rename( "Region" =Var1,  "Annotation" =Var2,"Population"=Freq)
	PopulationInRegion$Annotation <- factor(  PopulationInRegion$Annotation , levels = names(col_names))
	p3 <- PopulationInRegion %>% ggplot( aes( x= Region ,y = Population , fill= Annotation )) + geom_bar( stat = "identity", position = "fill") +
				 scale_fill_manual( values = col_names) + labs(title = "Population Plot" , subtitle = paste0(Plot_annotation) ) +
				 theme(axis.text.x = element_text(angle = 45, hjust = 1))  

	filePath <- paste0( "./PopulationPlot_",Plot_annotation,"_",name,".pdf")
    ggsave(file = filePath, plot = p3 , dpi=dpi, width=3, height=6)

    plot_resSplit  %>% as.data.frame %>% select(-c( total , adj_total , rate , adj_rate ) ) %>% write.table(.,file=paste0("MetadataAddRegion_",Plot_annotation,"_",name,".csv" ), quote=F, col.names=NA , sep="\t" )
    return(plot_resSplit)
}





