# Visium analysis

## Requirements
R (>=4.0; tested with 4.0.2 and 4.2.1)

The following R packages are required.
1. Seurat (v4)
2. hdf5r (for loding HDF5 files)

## Analysis
### Basic analysis of Visium data  
- [Visium_Seurat.R](./Visium_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat (Seurat 4.0.0 for IA, Seurat 4.3.0 for AIS/MIA).

### Integration of the data from serial sections 
