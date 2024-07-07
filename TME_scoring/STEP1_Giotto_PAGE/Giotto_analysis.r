library(Giotto)
library(patchwork)
library(ggplot2)
set.seed(1234)

SRanger <-"../demo_data/Spacerager/" # commandArgs(trailingOnly=TRUE)[1] 
Dir     <-"test" # commandArgs(trailingOnly=TRUE)[2] 
Pname   <-"test" # commandArgs(trailingOnly=TRUE)[3] 

#  5250    1000    3950    1000

temp_dir = paste0(Dir,"/fig")
dir.create(Dir)
dir.create(temp_dir) 

cellMaker <- read.table("../demo_data/PAGE_genesets.txt" , sep="\t" , header= T )

typeCol <- read.table("../demo_data/CelltypeCol.txt" , sep="\t" , header= T )
col_names <-  typeCol$col
names(col_names) <-  typeCol$CellType

instrs = createGiottoInstructions(save_dir = temp_dir , save_plot = TRUE , show_plot = FALSE)
GiottoOBJ = createGiottoVisiumObject(visium_dir = SRanger , expr_data = 'raw',
                                         png_name = 'tissue_lowres_image.png',
                                         gene_column_index = 2, instructions = instrs)

## subset on spots that were covered by tissue
metadata = pDataDT(GiottoOBJ)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
GiottoOBJ = subsetGiotto(GiottoOBJ, cell_ids = in_tissue_barcodes)

## filter
GiottoOBJ <- filterGiotto(gobject = GiottoOBJ,
                              expression_threshold = 1,
                              gene_det_in_min_cells = 1,
                              min_det_genes_per_cell = 10,
                              expression_values = c('raw'),
                              verbose = T)

## normalize
GiottoOBJ <- normalizeGiotto(gobject = GiottoOBJ, scalefactor = 10000, verbose = T)
## add gene & cell statistics

gene_metadata = fDataDT(GiottoOBJ)
cellMaker[ !cellMaker$Gene %in% gene_metadata$gene_ID  , ]
cellMaker <- cellMaker[ cellMaker$Gene %in% gene_metadata$gene_ID  , ]

sign_lsit <- list()
for(i in unique(cellMaker$CellType) ) { sign_lsit <- c( sign_lsit , list(cellMaker[cellMaker$CellType == i ,1 ] ))    }

signature_matrix = makeSignMatrixPAGE(sign_names = unique(cellMaker$CellType), sign_list =sign_lsit )
GiottoOBJ = runPAGEEnrich(gobject = GiottoOBJ, sign_matrix = signature_matrix , min_overlap_genes =2)
cell_types = colnames(signature_matrix)
spatCellPlot(gobject = GiottoOBJ,  spat_enr_names = 'PAGE',
             cell_annotation_values = cell_types,
             cow_n_col = 3,coord_fix_ratio = NULL, point_size = 0.75,
             save_param = list(save_name="7_b_spatcellplot_1", base_width = 12, base_height = 12) )

Annot <- as.data.frame(GiottoOBJ@spatial_enrichment$PAGE)
Annot$CellType <- apply( Annot,1,function(x){  names(x)[which.max(x)] })
GiottoOBJ@cell_metadata$CellType <-  Annot$CellType

p <- spatPlot2D(gobject = GiottoOBJ, cell_color  ="CellType"    , show_image = F, point_alpha = 1 , point_size = 2.3 ,
				axis_text =1 ,axis_title = 1, legend_symbol_size = 3, legend_text = 12, cell_color_code=col_names )
filePath <- paste0(Dir,"/Image_Annotation_", Pname ,".pdf")
ggsave(file = filePath, plot = p, dpi=100, width=5, height=5)

File <-paste0(Pname,"_PAGEscore.txt")   #"FFPE_No3B_PAGEscore.txt"
write.table( Annot   , file=File , sep="\t" , col.names=NA , row.names=T , quote=F)

