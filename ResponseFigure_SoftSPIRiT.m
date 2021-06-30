%-------------------------------------------------------------------------%
% Code for reproducing *** Figure S1 in Supporting Material *** for 
%     "Aliasing-free Reduced Field-of-View Parallel Imaging"
% Note: Run "CSPI_RedFOV_simulation" first for data preparation 
%-------------------------------------------------------------------------%
clc;clear;close all
addpath(strcat(pwd,'\data'))
addpath(genpath(strcat(pwd,'\utilities')))

% =========================================================================
% Data prepation
% =========================================================================
[Data_under_redFOV,Data_under_fullFOV,kCalib_fullFOV,kCalib_redFOV,Data_full] = ...
    CSPI_RedFOV_simulation;

% (reduced FOV) matrix size
[Ny_rFOV,Nz_rFOV,~] = size(Data_under_redFOV);

% (full FOV) matrix size
[Ny_fFOV,Nz_fFOV,~] = size(Data_under_fullFOV);


% =========================================================================
% Coil Compression
% =========================================================================
fprintf('Coil Compression ... \n');
Nc = 18;
sccmtx = calcSCCMtx(kCalib_fullFOV);
sccmtx_cc = sccmtx(:,1:Nc);

kCalib_redFOV = CC(kCalib_redFOV,sccmtx_cc);
kCalib_fullFOV = CC(kCalib_fullFOV,sccmtx_cc);

Data_under_redFOV = CC(Data_under_redFOV,sccmtx_cc);
Data_under_fullFOV = CC(Data_under_fullFOV,sccmtx_cc);
Data_full = CC(Data_full,sccmtx_cc);

% =========================================================================
% SPIRiT Kernel Calibration    
% =========================================================================
kSize = [7,7];      % SPIRiT kernel size
CalibTyk = 0.01;    % Tykhonov regularization in the calibration

disp('performing reduced FOV kernel calibration for SPIRiT');
tic;
kernel_redFOV = calibSPIRiT(kCalib_redFOV,kSize,Nc,CalibTyk);
GOP_rFOV = SPIRiT(kernel_redFOV,'image',[Ny_rFOV,Nz_rFOV]);
toc;


disp('performing full FOV kernel calibration for SPIRiT');
tic;
kernel_fullFOV = calibSPIRiT(kCalib_fullFOV,kSize,Nc,CalibTyk);
GOP_fFOV = SPIRiT(kernel_fullFOV, 'image',[Ny_fFOV,Nz_fFOV]);
toc;


% =========================================================================
% SPIRiT Conjugate Gradient Reconstruction    
% =========================================================================
nIterCG = 60;       % number of iteration; 

disp('CG-SPIRiT recon for reduced FOV acq + full FOV calibration & recon');
x0 = zeros(1,Ny_fFOV,Nz_fFOV,Nc);
y = reshape(Data_under_redFOV,[1,Ny_rFOV,Nz_rFOV,Nc]);
lambda = 1;
tic
[res_red2full_iterative,LSVEC] = cgSPIRiT_cartesian_img_v6(y, GOP_fFOV, nIterCG, lambda, x0);
toc
im_red2full_cartesian = sos(res_red2full_iterative);
as(squeeze(im_red2full_cartesian))


disp('CG-SPIRiT recon for full FOV acq & calibration & recon');
x0 = zeros(1,Ny_fFOV,Nz_fFOV,Nc);
y = reshape(Data_under_fullFOV,[1,Ny_fFOV,Nz_fFOV,Nc]);
lambda = 1;
tic
[res_fullFOV_cartesian] = cgSPIRiT_cartesian_img_v6(y, GOP_fFOV, nIterCG, lambda, x0);
toc
im_fullFOV_cartesian = sos(res_fullFOV_cartesian);
as(squeeze(im_fullFOV_cartesian))


disp('CG-SPIRiT recon for reduced FOV acq & calibration & recon');
x0 = zeros(1,Ny_rFOV,Nz_rFOV,Nc);
y = reshape(Data_under_redFOV,[1,Ny_rFOV,Nz_rFOV,Nc]);
lambda = 1;
tic
[res_redFOV_cartesian] = cgSPIRiT_cartesian_img_v6(y, GOP_rFOV, nIterCG, lambda, x0);
toc
im_redFOV_cartesian = sos(res_redFOV_cartesian);
as(squeeze(im_redFOV_cartesian))


% =========================================================================
% Post Soft-SENSE resolving the residual reduced FOV aliasing
% =========================================================================
reconksp_redFOV = fft2c(squeeze(res_redFOV_cartesian));

% A "Heavy" Coil Compression to only 4-6 channels because only 2-fold PI reconstruction
Nc = 6;   
sccmtx = calcSCCMtx(kCalib_fullFOV);
sccmtx_cc = sccmtx(:,1:Nc);
kCalib_fullFOV_4c = CC(kCalib_fullFOV,sccmtx_cc);
Data_redFOV = CC(reconksp_redFOV,sccmtx_cc);

% Full FOV CSM estimation of 4 channels
kSize = [5,5];  eigThresh_1 = 0.02;  eigThresh_2 = 0.95;
imSize = [Ny_fFOV,Nz_fFOV];
ecalib_map = espirit_calib2d(kCalib_fullFOV_4c,imSize,kSize,eigThresh_1,eigThresh_2,1);


% Mutiple-set CSMs used by Soft-SENSE
% nominal reduced FOV acceleration factors
Rfy = ceil(Ny_fFOV / Ny_rFOV);
Rfz = ceil(Nz_fFOV / Nz_rFOV);

% virtual FOV: multiple of reduced FOV and >= full FOV
Ny_vFOV = Ny_rFOV * Rfy;
Nz_vFOV = Nz_rFOV * Rfz;

% csm: full FOV zpad to virtual FOV
sens_map_vFOV = zpad(ecalib_map,[Ny_vFOV,Nz_vFOV,Nc]);

% shift
shift_amount_y = floor(Ny_vFOV/2) - floor(Ny_rFOV/2);
shift_amount_z = floor(Nz_vFOV/2) - floor(Nz_rFOV/2);
sens_map_shift = circshift(sens_map_vFOV,[-shift_amount_y,-shift_amount_z,0]);

% create (Rfy * Rfz) sens maps according to the fold-over process
sens_map_shift = reshape(sens_map_shift,[Ny_rFOV,Rfy,Nz_rFOV,Rfz,Nc]);
CSM_4sets = permute(sens_map_shift,[1,3,5,2,4]);


% Soft-SENSE to solve the aliasing due to reduced FOV only (Ry = 1 & Rz = 1)
Img0 = zeros(Ny_fFOV,Nz_fFOV);
[Img_red2full_twosteps,Gmap_red2full] = ...
    SoftSENSE_direct(Data_redFOV,CSM_4sets,Img0,1,1);

as(Img_red2full_twosteps)


% im_ref1 = sum(conj(squeeze(ecalib_map_fullFOV)).*ifft2c(Data_full),3);
im_ref2 = sos(ifft2c(Data_full));
diff_r2f = mat2gray(abs(im_ref2)) - mat2gray(abs(squeeze(im_red2full_cartesian)));
diff_full = mat2gray(abs(im_ref2)) - mat2gray(abs(squeeze(im_fullFOV_cartesian)));
diff_r2f_twosteps = mat2gray(abs(im_ref2)) - mat2gray(abs(Img_red2full_twosteps));

as(cat(2,diff_r2f,diff_full,diff_r2f_twosteps))
