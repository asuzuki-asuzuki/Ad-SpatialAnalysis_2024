# PhenoCycler analysis

## Requirements
1. QuPath (tested with 0.3.2)
2. QuPath StarDist extension (tested with 0.3.2)
3. R (>=4.0; tested with 4.3.3)

The following R packages are required.
1. Seurat (v5.1.0)

## Analysis
### Cell segmentation
1. A QPTIFF file (output file of PhenoCycler) is open by QuPath (version 0.3.2).
2. Regions are selected for tissue annotation.
3. Cell segmentation is performed using StarDist (QuPath StarDist extension, version 0.3.2) based on the DAPI signal (channel 0). The pixel size is set to 0.37 μm and the size of nucleus expansion is set to 5 µm. The script and model (stardist_cell_seg_model.pb) are provided by Akoya Biosciences.
4. After cell segmentation, mean pixel intensities for each marker (Cell: Mean) are exported as signal intensities in the CSV format.

### Analysis of segmented cells  
- [PhenoCycler_Seurat.R](./PhenoCycler_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat.

