# Visium analysis

## Requirements
R (>=4.0; tested with 4.0.2 and 4.2.1)  
Python (3.7 or later)  

The following R packages are required.
1. Seurat (v4)
2. hdf5r (for loding HDF5 files)
  
It generally takes only several tens of minutes for installlation of these packages.  
For conducting analyses, it usually takes short time (up to several hours). Memory requirements depend on the data size.  

## Analysis
### Basic analysis of Visium data  
- [Visium_Seurat.R](./Visium_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat.

### Visualization
- [Visium_Seurat_Violin.R](./Visium_Seurat_Violin.R): drawing violin plots of marker genes for each cluster by Seurat.
- [Visium_Seurat_ROI.R](./Visium_Seurat_ROI.R): drawing spatial plots of marker genes for each ROI by Seurat.
- [ModuleScore_Vizualization.R](./ModuleScore_Vizualization.R): calculating module scores using siganture genes and drawing spatial plots of the scores by Seurat.

### Merging the data of multiple serial sections 
- [FFPE_Merge_Seurat.R](./FFPE_Merge_Seurat.R): merging multiple Visium data by Seurat.

### Analysis of ligand-receptor interaction
- [COMMOT_LUAD3B_LigandReceptor_Calculation-checkpoint.ipynb](./COMMOT_LUAD3B_LigandReceptor_Calculation-checkpoint.ipynb): ligand-receptor communication analysis by COMMOT.
