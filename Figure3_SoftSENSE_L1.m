%-------------------------------------------------------------------------%
% Code for reproducing *** Figure 3 *** for 
%     "Aliasing-free Reduced Field-of-View Parallel Imaging"
% Note: Run "PI_RedFOV_simulation" first for data preparation 
%-------------------------------------------------------------------------%
clc;clear;close all

do_PICS_using_BART = 0;   % 0 - PI; 1 - PICS 

% =========================================================================
% Data prepation
% =========================================================================
[Data_under_redFOV,Data_under_fullFOV,kCalib_fullFOV,kCalib_redFOV,Data_full] = ...
    CSPI_RedFOV_simulation;

 
% (reduced FOV) matrix size
[Ny_rFOV,Nz_rFOV,~] = size(Data_under_redFOV);

% (full FOV) matrix size
[Ny_fFOV,Nz_fFOV,~] = size(Data_under_fullFOV);

% nominal reduced FOV acceleration factors
% i.e. times of fold-over caused by reduced FOV
Rfy = ceil(Ny_fFOV / Ny_rFOV);
Rfz = ceil(Nz_fFOV / Nz_rFOV);

% virtual FOV: multiple of reduced FOV and >= full FOV
Ny_vFOV = Ny_rFOV * Rfy;
Nz_vFOV = Nz_rFOV * Rfz;


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
% ESPIRiT Calibration of Coil Sensitivity Maps
% =========================================================================

if do_PICS_using_BART
    disp('ESPIRiT calibration of single-set full FOV CSM');
    kCalib_bart = reshape(zpad(kCalib_fullFOV,[Ny_fFOV,Nz_fFOV,Nc]),[1,Ny_fFOV,Nz_fFOV,Nc]);
    ecalib_map_fullFOV = bart('ecalib -m1 -k6 -c0.9 -S',kCalib_bart );
    load('/home/amax/ReducedFOVCode/data/reducedFOV/40ch/slice3/mask_roi.mat');
    ecalib_map_fullFOV = (ecalib_map_fullFOV) .* repmat(reshape(mask_roi,[1,Ny_fFOV,Nz_fFOV]),[1,1,1,Nc]);
    
    disp('ESPIRiT calibration of multiple-set reduced FOV CSM');
    kCalib_bart = reshape(zpad(kCalib_redFOV,[Ny_rFOV,Nz_rFOV,Nc]),[1,Ny_rFOV,Nz_rFOV,Nc]);
    ecalib_map_redFOV = bart('ecalib -m4 -k6 -c0.8 -S',kCalib_bart );
    
else
    
    disp('ESPIRiT calibration of single-set full FOV CSM');
    kSize = [6,6]; eigThresh_1 = 0.01; eigThresh_2 = 0.92;
    Nmaps = 1;
    ecalib_map = espirit_calib2d(kCalib_fullFOV,[Ny_fFOV,Nz_fFOV],kSize,eigThresh_1,eigThresh_2,Nmaps);
    ecalib_map_fullFOV = reshape(ecalib_map,[1,Ny_fFOV,Nz_fFOV,Nc]);
    load('/home/amax/ReducedFOVCode/data/reducedFOV/40ch/slice3/mask_roi.mat');
    ecalib_map_fullFOV = (ecalib_map_fullFOV) .* repmat(reshape(mask_roi,[1,Ny_fFOV,Nz_fFOV]),[1,1,1,Nc]);
    
    disp('ESPIRiT calibration of multiple-set reduced FOV CSM');
    kSize = [6,6];  eigThresh_1 = 0.01; eigThresh_2 = 0.88;
    Nmaps = 4;
    ecalib_map = espirit_calib2d(kCalib_redFOV,[Ny_rFOV,Nz_rFOV],kSize,eigThresh_1,eigThresh_2,Nmaps);
    ecalib_map_redFOV = reshape(ecalib_map,[1,Ny_rFOV,Nz_rFOV,Nc,Nmaps]);
end

% =========================================================================
% Soft-SENSE reconstructions
% =========================================================================
nIterCG = 30;
Nx = 1;
L1w = 0.01;

disp('Soft-SENSE (4-set CSM) recon for 2x2 PI + 2.25x reduced FOV acq + full FOV calib & recon');
CSM_fullFOV_zpad = zpad(ecalib_map_fullFOV, [Nx,Ny_rFOV*2,Nz_rFOV*2,Nc]);

CSM_4sets = permute(reshape(circshift(CSM_fullFOV_zpad, [0,Ny_rFOV/2,0,0]),...
    [Nx,Ny_rFOV,2,Nz_rFOV*2,Nc]), [1,2,4,5,3]);

CSM_4sets = permute(reshape(circshift(CSM_4sets, [0,0,Nz_rFOV/2,0,Nz_rFOV/2]), ...
    [Nx,Ny_rFOV,Nz_rFOV,2,Nc,2]), [1,2,3,5,6,4]);

CSM_4sets = reshape(CSM_4sets,[Nx,Ny_rFOV,Nz_rFOV,Nc,4]);

x0 = zeros(1,Ny_rFOV,Nz_rFOV,1,4);
ksp_red = reshape(Data_under_redFOV,[1,Ny_rFOV,Nz_rFOV,Nc]);
mask = abs(ksp_red) > 0;
tic
if do_PICS_using_BART
   bart_pics_command = sprintf('pics -l1 -i%d -r%f -m -S',nIterCG,L1w);
   im_red2full_4sets = bart(bart_pics_command,ksp_red,CSM_4sets);
else
   im_red2full_4sets = cgSoftSENSE_cartesian_v3(ksp_red, x0, CSM_4sets,mask, nIterCG);
end
toc


img_fullFOV_shift = reshape(im_red2full_4sets,[Nx,Ny_rFOV,Nz_rFOV,2,2]);

img_fullFOV_shift = permute(img_fullFOV_shift,[1,2,4,3,5]);

img_fullFOV_shift = flip(img_fullFOV_shift,5);
tmp = circshift(img_fullFOV_shift(:,:,:,:), [0,0,0,Nz_rFOV/2]);
img_fullFOV_shift = permute(tmp,[1,4,2,3]);
tmp = circshift(img_fullFOV_shift(:,:,:), [0,0,Ny_rFOV/2]);
img_final = permute(tmp,[1,3,2]);

im_red2full = squeeze(crop(img_final,[Nx,Ny_fFOV,Nz_fFOV]));
as(im_red2full)

disp('Soft-SENSE (2-set CSM) recon for 2x2 PI + 2.25x reduced FOV acq + full FOV calib & recon');
CSM_2sets = zeros(Nx,Ny_rFOV,Nz_rFOV,Nc,2);
CSM_2sets(:,:,:,:,1) = CSM_4sets(:,:,:,:,3);
CSM_2sets(:,:,:,:,2) = CSM_4sets(:,:,:,:,1) + CSM_4sets(:,:,:,:,4);

x0 = zeros(1,Ny_rFOV,Nz_rFOV,1,2);
tic
if do_PICS_using_BART
   bart_pics_command = sprintf('pics -l1 -i%d -r%f -m -S',nIterCG,L1w);
   im_red2full_2sets = bart(bart_pics_command,ksp_red,CSM_2sets);
else
   im_red2full_2sets = cgSoftSENSE_cartesian_v3(ksp_red, x0, CSM_2sets,mask, nIterCG);
end
toc
as(squeeze(im_red2full_2sets))


disp('Soft-SENSE recon for 2x2 PI + 2.25x reduced FOV acq & calib & recon');
x0 = zeros(1,Ny_rFOV,Nz_rFOV,1,size(ecalib_map_redFOV,5));
ksp = reshape(Data_under_redFOV,[1,Ny_rFOV,Nz_rFOV,Nc]);
mask = abs(ksp) > 0;

tic
if do_PICS_using_BART
    bart_pics_command = sprintf('pics -l1 -i%d -r%f -m -S',nIterCG,L1w);
    im_redFOV = bart(bart_pics_command,ksp_red,ecalib_map_redFOV);
else
    im_redFOV = cgSoftSENSE_cartesian_v3(ksp, x0, ecalib_map_redFOV, mask, nIterCG);
end
toc
as(squeeze(im_redFOV))


disp('Soft-SENSE recon for 3x3 full FOV acq & calib & recon');
x0 = zeros(1,Ny_fFOV,Nz_fFOV,1,1);
ksp = reshape(Data_under_fullFOV,[1,Ny_fFOV,Nz_fFOV,Nc]);
mask = abs(ksp) > 0;

tic
if do_PICS_using_BART
    bart_pics_command = sprintf('pics -l1 -i%d -r%f -m -S',nIterCG,L1w);
    im_fullFOV = bart(bart_pics_command,ksp,ecalib_map_fullFOV);
else
    im_fullFOV = cgSoftSENSE_cartesian_v3(ksp, x0, ecalib_map_fullFOV, mask, nIterCG);
end
toc
as(squeeze(im_fullFOV))


% im_ref = sum(conj(squeeze(ecalib_map_fullFOV)).*ifft2c(Data_full),3);
% diff_r2f = mat2gray(abs(im_ref)) - mat2gray(abs(im_red2full));
% diff_full = mat2gray(abs(im_ref)) - mat2gray(abs(squeeze(im_fullFOV)));
% 
% as(cat(2,diff_r2f,diff_full))

