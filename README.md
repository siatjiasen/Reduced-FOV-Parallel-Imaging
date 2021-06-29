# Reduced-FOV-Parallel-Imaging

Code for reproducing the Figures in 
"Aliasing-free Reduced Fielf-of-View Parallel Imaging"

This package relies on the two external packages:

1. arrayShow - displaying multi-dimensional image array (https://github.com/tsumpf/arrShow)
2. ESPIRiT   - MATLAB ESPIRiT calibration and utility functions (https://people.eecs.berkeley.edu/~mlustig/software/SPIRiT_v0.3.tar.gz)
3. BART      - C ESPIRiT calibration and Compressed Sensing Parallel Imaging reconstruction (https://github.com/mrirecon/bart)

Figure 1 - 2D imaging with 1D Reduced FOV Parallel Imaging

Figure 2 - 3D imaging with 2D Reduced FOV Parallel Imaging 

Figure 3 - 3D imaging with 2D Reduced FOV Compressed Sensing

Figure 4 - 2D imaging with 1D Reduced FOV Wave Encoding

Figure 5 - 3D imaging with 2D Reduced FOV Wave Encoding

SoftSENSE_Direct      - SoftSENSE reconstruction for regular downsampling by pixel-by-pixel direct inversion

cgSoftSENSE_Cartesian - SoftSENSE reconstruction for 2D/3D Cartesian imaging with/without reduced FOV imaging

cgSoftSENSE_Wave      - SoftSENSE reconstruction for 2D/3D Wave imaging with/without reduced FOV imaging

cgSPIRiT_Cartesian    - SPIRiT    reconstruction for 2D/3D Cartesian imaging with/without reduced FOV imaging (SPIRiT image model)

cgSPIRiT_Wave         - SPIRiT    reconstruction for 2D/3D Wave imaging with/without reduced FOV imaging

