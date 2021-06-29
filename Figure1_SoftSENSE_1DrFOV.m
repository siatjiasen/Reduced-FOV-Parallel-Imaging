%-------------------------------------------------------------------------%
% Code for reproducing *** Figure 2 *** for 
%     "Aliasing-free Reduced Field-of-View Parallel Imaging"
% Note: Run "PI_RedFOV_simulation" first for data preparation 
%-------------------------------------------------------------------------%
clc;clear;close all

addpath(strcat(pwd,'\data'))
addpath(genpath(strcat(pwd,'\utilities')))

using_bart_csm = 0;   % 0 - using ESPIRiT implemented in MATLAB, slow; 

%  reduced FOV imaging acceleration
%  Net acceleration factor = 2.25

Ry_redFOV = 1;
Rz_redFOV = 1.5;


% parallel imaging acceleration factor for reduced FOV imaging
% Net acceleration factor :  3 = 2x1.5
Ry_PI_redFOV = 1;
Rz_PI_redFOV = 2;

% parallel imaging acceleration factor for full FOV imaging
% Net acceleration factor : 3
Ry_PI_fullFOV = 1;
Rz_PI_fullFOV = 3;


% =========================================================================
% Data prepation
% =========================================================================
[Data_under_redFOV,Data_under_fullFOV,kCalib_fullFOV,kCalib_redFOV] = ...
    PI_RedFOV_simulation(Ry_redFOV,Rz_redFOV,Ry_PI_redFOV,Rz_PI_redFOV,Ry_PI_fullFOV,Rz_PI_fullFOV);

% (reduced FOV) kspace matrix size
[Ny_rFOV,Nz_rFOV,~] = size(Data_under_redFOV);

% (full FOV) kspace matrix size
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


% =========================================================================
% ESPIRiT Calibration of Coil Sensitivity Maps 
%       Full FOV calibration     : single set      (-m1)
%       Reduced FOV calibraation : multiple sets   (-m4)
% =========================================================================
fprintf('ESPIRiT Calibration ... \n');
if using_bart_csm
    kCalib_bart = reshape(zpad(kCalib_fullFOV,[Ny_fFOV,Nz_fFOV,Nc]),[1,Ny_fFOV,Nz_fFOV,Nc]);
    ecalib_map_fullFOV = bart('ecalib -a -m1 -k6',kCalib_bart);
    load('./data/mask_roi.mat');
    CSM_fullFOV = squeeze(ecalib_map_fullFOV) .* repmat(mask_roi,[1,1,Nc]);
    
    kCalib_bart = reshape(zpad(kCalib_redFOV,[Ny_rFOV,Nz_rFOV,Nc]),[1,Ny_rFOV,Nz_rFOV,Nc]);
    if Ry_redFOV == 1 || Rz_redFOV == 1
        ecalib_map_redFOV  = bart('ecalib -a -m2  -k6',kCalib_bart);
    else
        ecalib_map_redFOV = bart('ecalib -a -m4 -k6',kCalib_bart);
    end
    CSM_ESPIRiT = squeeze(ecalib_map_redFOV);
    
else
    
    kSize = [6,6]; eigThresh_1 = 0.01; eigThresh_2 = 0.92;
    Nmaps = 1;
    ecalib_map = espirit_calib2d(kCalib_fullFOV,[Ny_fFOV,Nz_fFOV],kSize,eigThresh_1,eigThresh_2,Nmaps);
    ecalib_map_fullFOV = reshape(ecalib_map,[1,Ny_fFOV,Nz_fFOV,Nc]);
    load('./data/mask_roi.mat');
    CSM_fullFOV = squeeze(ecalib_map_fullFOV) .* repmat(mask_roi,[1,1,Nc]);
    
    kSize = [6,6];  eigThresh_1 = 0.01; eigThresh_2 = 0.88;
    if Ry_redFOV == 1 || Rz_redFOV == 1
        Nmaps = 2;
    else
        Nmaps = 4;
    end
    ecalib_map = espirit_calib2d(kCalib_redFOV,[Ny_rFOV,Nz_rFOV],kSize,eigThresh_1,eigThresh_2,Nmaps);
    CSM_ESPIRiT = squeeze(ecalib_map(:,:,:,end:-1:1));
    
end

% =========================================================================
% Create mulitple-set CSM for soft-SENSE according to the fold-over of
% reduced FOV imaging
% =========================================================================
% zpad to virtual FOV
sens_map_vFOV = zpad(CSM_fullFOV,[Ny_vFOV,Nz_vFOV,Nc]);

% FOV shifting
shift_amount_y = floor(Ny_vFOV/2) - floor(Ny_rFOV/2);
shift_amount_z = floor(Nz_vFOV/2) - floor(Nz_rFOV/2);
sens_map_shift = circshift(sens_map_vFOV,[-shift_amount_y,-shift_amount_z,0]);

% fold-over due to reduced FOV imaging
sens_map_shift = reshape(sens_map_shift,[Ny_rFOV,Rfy,Nz_rFOV,Rfz,Nc]);
sens_map_shift = permute(sens_map_shift,[1,3,5,2,4]);
CSM_SoftSENSE = sens_map_shift;


% =========================================================================
% Soft-SENSE Recon (Direct Inversion, only works for regular downsampling)
% =========================================================================
Img0 = zeros(Ny_fFOV,Nz_fFOV);
[Img_red2full,Gmap_red2full] = ...
    SoftSENSE_direct(Data_under_redFOV,CSM_SoftSENSE,Img0,Ry_PI_redFOV,Rz_PI_redFOV);

Img0 = zeros(Ny_fFOV,Nz_fFOV);
[Img_fullFOV,Gmap_fullFOV] = ...
    SoftSENSE_direct(Data_under_fullFOV,CSM_fullFOV,Img0,Ry_PI_fullFOV,Rz_PI_fullFOV);

Img0 = zeros(Ny_fFOV,Nz_fFOV);
[Img_redFOV,Gmap_redFOV] = ...
    SoftSENSE_direct(Data_under_redFOV,CSM_ESPIRiT,Img0, Ry_PI_redFOV,Rz_PI_redFOV);
Img_redFOV = zpad(crop(Img_redFOV,[Ny_rFOV,Nz_rFOV]),[Ny_fFOV,Nz_fFOV]);
Gmap_redFOV = zpad(crop(Gmap_redFOV,[Ny_rFOV,Nz_rFOV]),[Ny_fFOV,Nz_fFOV]);

Img0 = zeros(Ny_rFOV,Nz_rFOV);
[Img_redFOV_1map,Gmap_redFOV_1map] = ...
    SoftSENSE_direct(Data_under_redFOV,CSM_ESPIRiT(:,:,:,1),Img0, Ry_PI_redFOV,Rz_PI_redFOV);
Img_redFOV_1map = zpad(crop(Img_redFOV_1map,[Ny_rFOV,Nz_rFOV]),[Ny_fFOV,Nz_fFOV]);
Gmap_redFOV_1map = zpad(crop(Gmap_redFOV_1map,[Ny_rFOV,Nz_rFOV]),[Ny_fFOV,Nz_fFOV]);

% =========================================================================
% Results Display
% =========================================================================
g_max_r2f = max(Gmap_red2full(:));
g_max_full = max(Gmap_fullFOV(:));
g_max_red = max(Gmap_redFOV(:));
figure;
imshow(cat(2,mat2gray(abs(Img_red2full)),mat2gray(abs(Img_fullFOV)),mat2gray(abs(Img_redFOV))));
title('SoftSENSE:PI 2x2+ReducedFOV        SENSE:PI 3x3+Full FOV          ESPIRiT:PI 2x2+ReducedFOV ')

figure;
g_disp_max = max([g_max_r2f,g_max_full,g_max_red]);
imshow(cat(2,Gmap_red2full,Gmap_fullFOV,Gmap_redFOV),[0,3.5]);
colormap(jet(256)); colorbar;
title({'SoftSENSE:PI 2x2+ReducedFOV        SENSE:PI 3x3+Full FOV          ESPIRiT:PI 2x2+ReducedFOV ';...
    strcat(strcat('Gmax:',num2str(g_max_r2f)), ...
    strcat('                        Gmax:',num2str(g_max_full)),...
    strcat('                                 Gmax:',num2str(g_max_red)))});



