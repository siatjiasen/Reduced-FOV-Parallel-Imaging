# Reduced-FOV-Parallel-Imaging

Code for reproducing the Figures in 
"Aliasing-free Reduced Fielf-of-View Parallel Imaging"

This package relies on the following external packages:

1. arrayShow - displaying multi-dimensional image array (https://github.com/tsumpf/arrShow)
2. ESPIRiT   - MATLAB ESPIRiT calibration and utility functions (https://people.eecs.berkeley.edu/~mlustig/software/SPIRiT_v0.3.tar.gz)
3. BART      - C ESPIRiT calibration and Compressed Sensing Parallel Imaging reconstruction (https://github.com/mrirecon/bart)

Figure 1 - 2D imaging with 1D Reduced FOV Parallel Imaging (MATALB)

Figure 2 - 3D imaging with 2D Reduced FOV Parallel Imaging (MATLAB, a 2D slice extracted from 3D kspace)

Figure 3 - 3D imaging with 2D Reduced FOV Compressed Sensing (MATLAB, a 2D slice extracted from 3D kspace)

Figure 4 - 2D imaging with 1D Reduced FOV Wave Encoding (MATLAB)

Figure 5 - 3D imaging with 2D Reduced FOV Wave Encoding (BART)  
           

Data Description: 
All experimental datasets (Total 3.4 GB) can be downloaded from https://share.weiyun.com/LoDm1LAe . 
The data should be stored into the folder named as /data.

The experimental dataset for Figure 1,2,3 can be downloaded separately from https://share.weiyun.com/qWy1E9sC (72 MB).

The experimental dataset for Figure 4 can be downloaded separately from https://share.weiyun.com/7rC5jbLI (150 MB).

The experimental dataset for Figure 5 including two 3D data of size of 2.2 GB (PI 2x2) and 1.2 GB (reduced FOV 2x + PI 2x) respectively.


The implementations of SoftSENSE:

SoftSENSE_Direct      - SoftSENSE reconstruction for regular downsampling by pixel-by-pixel direct inversion

cgSoftSENSE_Cartesian - SoftSENSE reconstruction for 2D/3D Cartesian imaging with/without reduced FOV imaging

cgSoftSENSE_Wave      - SoftSENSE reconstruction for 2D/3D Wave imaging with/without reduced FOV imaging

cgSPIRiT_Cartesian    - SPIRiT    reconstruction for 2D/3D Cartesian imaging with/without reduced FOV imaging

cgSPIRiT_Wave         - SPIRiT    reconstruction for 2D/3D Wave imaging with/without reduced FOV imaging

