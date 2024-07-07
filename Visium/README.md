# Visium analysis

## Requirements
R (>=4.0; tested with 4.0.2 and 4.2.1)

The following R packages are required.
1. Seurat (v4)
2. hdf5r (for loding HDF5 files)

## Analysis
### Basic analysis of Visium data  
- [Visium_Seurat.R](./Visium_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat.

### Visualization
- [Visium_Seurat_Violin.R](./Visium_Seurat_Violin.R): drawing violin plots of marker genes for each cluster by Seurat.
- [Visium_Seurat_ROI.R](./Visium_Seurat_ROI.R): drawing spatial plots of marker genes for each ROI by Seurat.

### Merging the data of multiple serial sections 
- [FFPE_Merge_Seurat.R](./FFPE_Merge_Seurat.R): merging multiple Visium data by Seurat.

### Analysis of ligand-receptor interaction
Under preparation
