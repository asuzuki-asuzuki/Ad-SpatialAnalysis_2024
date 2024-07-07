library(dplyr)
library(tidyverse)
library(ggplot2)
library(SPATA2)
library(gggibbous)

# Make_SCCR(spata_obj= spata_obj ,PAGE_inf=PAGE_inf ,project=project  ,trajectory_name=trajectory_name , slope_min=slope_min,region_length=region_length ,Annotation ) 
# # Parameters.
# - spata_obj ... SPATA object after calculating trajectories.
# - PAGE_inf  ... PAGE score matorix
# - project ... project name . this name used by ouput name
# - trajectory_name ... Trajectory name given as input to the SPTAT2 createTrajectoryManually() run.
# - slope_min == slope_min ... The lowest slope value in the distribution of trajectory scores calculated by SPTAT2 that is defined as having changed.
# - region_lengthPct ...  Defines the length of the region required to be defined as changed.
#                         Greater than the total length of the region multiplied by the value of "region_lengthPct". 
#                         (If region_length is defined, this value is ignored).
# - span ... loess span
# - overlap ... 
# - region_length ... 
# - Annotation ... Plot PAGE Category

# - pw ... width of ouput pdf
# - ph ... height of ouput pdf

Make_SCCR <- function ( spata_obj,
				PAGE_inf,
				project,
				trajectory_name,
				slope_min, 
				region_lengthPct ,
				span =0.2 ,
				overlap = 5 ,
				pw= 8,
				ph= 8,
				region_length = NULL,
				Annotation = NULL 
) {

	#%%%%%
	# -  Get Spatial Trajetory data , and plot by SPTAT2
	#%%%%%
	p1_lists <- list() # Trajectory Plot , SPTAT2 plotTrajectory() function.
	p2_lists <- list() # Heatmap of Each Anonotation PAGE score , nomlized by Trajectory. SPATA2 plotTrajectoryFeatures() function.
	tjk_mat <- c() 
	if(is.null(Annotation)) {
		Annotation = colnames( PAGE_inf)
	}

	numplot= length(Annotation)
	plot_ncol =4
	plot_nrow = ceiling( numplot/plot_ncol)

	for (i in c(1: (length(Annotation)  ) )  ) {
		Class =  as.character(Annotation[i])

		# plotTrajectory
		p1<- plotTrajectory(
					object = spata_obj, 
					trajectory_name = trajectory_name , 
					pt_alpha =1 ,
					pt_alpha2 =1, 
					sgmt_clr = "skyblue",
				    color_by = Class , 
					smooth = TRUE, 
					display_image = FALSE  
			 ) + 
			 legendTop() + 
			 theme(text = element_text(size = 12))

		# plotTrajectoryFeatures
		p2<- plotTrajectoryFeatures(
					object = spata_obj, 
					trajectory_name = trajectory_name , 
					features = Class ,
					smooth_method = "loess", 
					smooth_span = span,
					display_trajectory_parts = TRUE,
					 smooth_se = TRUE 
			 )  + 
			 theme(text = element_text(size = 16))
		# make matrix 
		if(i == 1 ){
			tjk_mat <-as.data.frame(p2$data)[,c(1,2,3,5)]
		} else {
			tjk_mat <- cbind( tjk_mat , as.data.frame(p2$data)[5] )
		}
		# append plot data
		p1_lists <- c( p1_lists , list(p1) )
		p2_lists <- c( p2_lists , list(p2) )
	}

	colnames(tjk_mat)[4:(ncol(tjk_mat))] <- as.vector(Annotation)
	# output , Matrix of each PAGE score at each trajectory.
	outputfile <-  paste0( "trajectory_score_matrix_", project,"_", trajectory_name ,".txt")
	write.table( tjk_mat , file=  outputfile , quote=F, col.names=NA , sep="\t")

	# save  Trajectory plot list
	filePath <- paste0( "./ALL_PAGE_Trajectory_", project,"_", trajectory_name ,".pdf")
    ggsave(file = filePath, plot = wrap_plots(p2_lists) , dpi=100, width=16, height=12)

	# save  heatmap  of PAGE score
	filePath <- paste0( "./ALL_PAGE_Heatmap_"   , project,"_", trajectory_name ,".pdf")
    ggsave(file = filePath, plot = wrap_plots(p1_lists) , dpi=100, width=16, height=12)

	#%%%%%
	# Searching for region of PAGE score change.
	#%%%%%

	#  if region_length , defins region length from region_lengthPct.
	if(is.null(region_length)){
		region_length = as.integer(nrow(tjk_mat) * region_lengthPct)
	}

	# Dectect SCCR
	pdf(paste0("./InflectionPoint_Trajectory_ScorePlot", project,"_", trajectory_name ,".pdf") , height= plot_nrow * 4 ,width = plot_ncol *4 )
	par(mfrow = c(plot_nrow, plot_ncol)) ###

	# list of  slope changed point.
    change_region_list <-  matrix( nrow = 0, ncol = 6)
    colnames( change_region_list ) = c("x_start","x_end","pos","Annotaion","sid","slope")

    # 
	for (i in c(1: length(Annotation))  ) {
		Class = as.character(Annotation[i])

		use_mat <- tjk_mat[,c("trajectory_order",Class) ]
		colnames(use_mat) <- c("trajectory_order","val")

		# calculation loess fitting.
		loessRes <-loess(formula =  val ~ trajectory_order, use_mat , model = TRUE,span = span ,method = "loess" )
		# get InflectionPoint list 
		smoothRes <- predict(loessRes)
		infl <- c(TRUE, diff(diff(smoothRes)>0)!=0,TRUE)
		infl <- c(TRUE, diff(diff(smoothRes)>0)!=0,TRUE)

		# plot row PAGE ccore on trajectory.
		plot(use_mat$trajectory_order, use_mat[,2], pch=19, main= paste0(Class) ,cex=0.7, xlab="", ylab="")
		# plot fitted line.
		lines(smoothRes, x= use_mat$trajectory_order, col='red')
		
		# Get Start point of each trajectory Part and plotting line 
		ch_trj <-  which(  tjk_mat$trajectory_part_order==1)
		for (n in c(1:length(ch_trj)) ){
			lines( y=c ( -1 , 1.5 ) ,x=c( ch_trj[n], ch_trj[n] ) , lty = "dotted", col = "grey" )
		}

		# plot InflectionPoint by yelow
		points(use_mat[infl, "trajectory_order" ], smoothRes[infl ], col="yellow",pch =20, cex =1)
		infl_point <- cbind( y= smoothRes[which( infl ,TRUE )],  x=which( infl ,TRUE ) )

		# Calculates the slope between two inflection points and extracts the area between the two points where the threshold is exceeded
		# and the width of the area is greater than the threshold.
		for (k in c(2:nrow(infl_point))) {
			va <-  ( infl_point[k,"y"] -infl_point[k-1,"y"]    ) /( infl_point[k,"x"] -infl_point[k-1,"x"]    ) # slope of 2 InflectionPoint
			x_1 <- infl_point[(k-1),"x"] 
			x_2 <- infl_point[k,"x"]        
			if( abs(va) > slope_min & (x_2 - x_1) > region_length   ){
				if( va > slope_min ) {
					polygon(x=c(x_1,x_1,x_2,x_2), y=c(0,1,1,0), col="#FF990022", border=F) 
					plus <- "+" # Positive Correlation with trajectory
				} else  {
					polygon(x=c(x_1,x_1,x_2,x_2), y=c(0,1,1,0), col="#0000FF22", border=F)
					plus <- "-" # Negative Correlation with trajectory
				}
				#sid =   sum(change_region_list[,"pos"] == plus &  change_region_list[,"Annotaion"] == Class) + 1
				sid =   sum( change_region_list[,"Annotaion"] == Class) + 1
				change_region_list <- rbind( change_region_list , c(x_1,x_2,plus,Class,sid,va  ) )
			}
		}
	}
	dev.off()

	#%%%%%
	# assotiation between celltypes.
	#%%%%%

	change_region_list <-  as.data.frame( change_region_list)
	change_region_list[,c(1,2,5,6)] <- data.frame( lapply(change_region_list[,c(1,2,5,6)] , as.numeric))
	print(change_region_list)

	ch_region_mat <-  matrix( nrow = 0, ncol = 15)
	colnames(ch_region_mat) <- c( colnames(change_region_list),"slope_ov" ,  paste0(colnames(change_region_list),"2"), "slope_ov2" , "overlap" )

	# Examine the relationship between two annotation categories
	for ( i in c(1: nrow(change_region_list)  ) ) {

		x_start <- change_region_list[i,"x_start"]
		x_end   <- change_region_list[i,"x_end"  ] 

		for (k in c(1:nrow(change_region_list))) {	
			if( change_region_list[i,"Annotaion"] ==  change_region_list[k,"Annotaion"] ){
				next
			}
			x2_start <- change_region_list[k,"x_start"]
			x2_end   <- change_region_list[k,"x_end"  ]

			# Check for overlap between two categories
			if ( !(  x_end <= x2_start |  x2_end <= x_start )  ){
				ov_x_start = ifelse( x2_start > x_start , x2_start  , x_start )
				ov_x_end   = ifelse( x2_end   < x_end   , x2_end    , x_end   )
				ov_lap     = (ov_x_end  - ov_x_start +1  ) /  nrow(tjk_mat) 

				# check ovelatp length
				if( ( ov_x_end - ov_x_start +1  ) >=  overlap  ){
					slop1 <- ( tjk_mat[ov_x_end,  change_region_list[i,"Annotaion"] ] - tjk_mat[ov_x_start, change_region_list[i,"Annotaion"] ]   ) / ( ov_x_end - ov_x_start  )
					slop2 <- ( tjk_mat[ov_x_end,  change_region_list[k,"Annotaion"] ] - tjk_mat[ov_x_start, change_region_list[k,"Annotaion"] ]   ) / ( ov_x_end - ov_x_start  )
					ch_region_mat <-rbind (ch_region_mat , c( change_region_list[i,1:6],slop1 , change_region_list[k,1:6],slop2,ov_lap  ))
				} else {
					ch_region_mat <-rbind (ch_region_mat , c( change_region_list[i,1:6], 0 , change_region_list[k,1:6], 0 , 0      ))
				}
			} else {
				ch_region_mat <-rbind (ch_region_mat , c( change_region_list[i,1:6], 0 , change_region_list[k,1:6], 0 , 0      ))
			}
		}
	}

	filename = paste0( "ScoreChangedRegoin_", project,"_", trajectory_name ,".txt")
	write.table(  ch_region_mat , file= filename , quote=F, col.names=NA , sep="\t")
	ch_region_mat <- ch_region_mat %>% as.data.frame()
	ch_region_mat [,c(1,2,5,6,7,8,9,12,13,14,15)] <- data.frame( lapply(ch_region_mat[,c(1,2,5,6,7,8,9,12,13,14,15)] , as.numeric))
	# x_start      x_end      pos     Annotaion       sid     slope   slope_ov  x_start2     x_end2     pos2    Annotaion2      sid2    slope2  slope_ov2 overlap

	# modified matrix
	forPlot_ch_mat <-ch_region_mat %>% 
		mutate( 
			NewID = paste(!!!rlang::syms(c("Annotaion","sid","pos")), sep=":") ,
			NewID2 = paste(!!!rlang::syms(c("Annotaion2","sid2","pos2")), sep=":") 
		) %>%
		select( c("NewID","NewID2","slope","slope2","overlap"  )) %>%
		pivot_longer( cols=c("slope","slope2"),names_to ="type", values_to ="slope" )  %>%
		mutate( 
			start = recode( type  , slope = -3.141593 ,slope2 = 0) , 
			end=recode( type  , slope = 0 , slope2 =  3.141593) , 
			right= recode( type  , slope = FALSE , slope2 =  TRUE)
		)

	# modified order
	defLevelsOfIDs <- unique(forPlot_ch_mat$NewID)
    splitedLevelsOfIDs <- as.data.frame(  t(matrix(unlist(strsplit(defLevelsOfIDs, ":")) , nrow=3 )) )
	splitedLevelsOfIDs$V1 <-  factor( splitedLevelsOfIDs$V1 , levels = Annotation) 
    sorted_splitedLevelsOfIDs <-with( splitedLevelsOfIDs , splitedLevelsOfIDs[order(V1, V3, V2),])
	NewOrder <- apply(sorted_splitedLevelsOfIDs,1,function(x){ return(paste(as.matrix(x), collapse = ":")) })
	forPlot_ch_mat$NewID  <- factor( forPlot_ch_mat$NewID  , levels = NewOrder)
	forPlot_ch_mat$NewID2 <- factor( forPlot_ch_mat$NewID2 , levels = NewOrder)
	
	forPlot_ch_mat <- forPlot_ch_mat %>%	
		mutate(
			x = match( NewID , levels(NewID)  ) ,
			y = match( NewID2 , levels(NewID))
		)

	# bubble plotting
	p<- ggplot(forPlot_ch_mat) +
		geom_moon( 
			aes(
				x = x,
				y = y,
				ratio = 0.5,
				right = right,
				size = overlap ,
				fill = slope
			),
			key_glyph = draw_key_full_moon
		) +
		scale_fill_gradient2(midpoint=0 , low="blue", mid="white",high="red" , limits =c(-0.08 , 0.08) , oob = scales::squish) +
		scale_size("size", range = c(1,10) ,   breaks =c(0,0.05,0.10,0.15,0.20,0.25,0.30) ,limit=c(0,0.3) ) +
	    scale_x_continuous(breaks= c(1: length(unique (forPlot_ch_mat$NewID))) ,labels=c(levels(forPlot_ch_mat$NewID))) +
	    scale_y_continuous(breaks= c(1: length(unique (forPlot_ch_mat$NewID))) ,labels=c(levels(forPlot_ch_mat$NewID))) +
	    xlab("") + ylab("") +
	    coord_fixed() +
	    theme_bw() +
	    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 , size = rel(1))) 

	filePath <- paste0("TJK_ovarlap_", project,"_", trajectory_name ,".pdf")
	ggsave(file = filePath, plot = p, dpi=100, width=pw, height=ph)
	return( list(change_region_list,tjk_mat) )
}


Fun_lotate <- function( loc , do){
    c = do *pi/180
    x = loc[1]
    y = loc[2]
    x2 = x*cos(c) - y * sin(c)
    y2 = x*sin(c) + y * cos(c)
    return(c(x2,y2))
}


