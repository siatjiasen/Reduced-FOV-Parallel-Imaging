%-------------------------------------------------------------------------%
% Code for reproducing *** Figure 5 *** for
%     "Aliasing-free Reduced Field-of-View Parallel Imaging"
%-------------------------------------------------------------------------%
clc;clear;close all
addpath(strcat(pwd,'/data'))
addpath(genpath(strcat(pwd,'/utilities')))

iscoil_compression = 1;
isremove_roos_acs = 1;

%--------------------------------------------------------------------------
% Step 0. rawdata (imaging & CSM + PSF calibration scan) reading
% -------------------------------------------------------------------------
%---- Full FOV calibration data

load ./data/WaveRefScan_fullFOV.mat

PsfYZ_fullFOV = wavepsf_calibration(Data_Py,Data_WavePy,Data_Pz,Data_WavePz);
ACS_fullFOV = squeeze(Data_ACS(:,:,:,:,1));

%----  Reduced FOV imaging data
RAW_redFOV = readcfl('./data/rawdata_redFOV'); % data is accelerated by rFOV 2x (PE1) + PI 2x (PE2)
% RAW_2x2 = readcfl('./data/rawdata_2x2');     % data is accelerated by PI 2x (PE1) + PI 2x (PE2)


acsdata = ACS_fullFOV;
rawdata = RAW_redFOV;

Nx_os = size(PsfYZ_fullFOV,1);
Ny_fullFOV = size(PsfYZ_fullFOV,2);
Nz_fullFOV = size(PsfYZ_fullFOV,3);

Ny_redFOV = size(rawdata,2);
Nz_redFOV = size(rawdata,3);

%--------------------------------------------------------------------------
% Step 1: ACS preparation
%--------------------------------------------------------------------------
nRO_os = size(acsdata,1);
if (isremove_roos_acs == 1)
    nRO = nRO_os/2;          % wave readout oversampling ratio
    acs_crop = fftshift(ifft(fftshift(acsdata,1),[],1),1);
    acs_crop = crop(acs_crop,[nRO,24,24,size(acsdata,4)]);
    Data_ACS = fftshift(fft(fftshift(acs_crop,1),[],1),1);
else
    Data_ACS = acsdata;
    nRO = nRO_os;
end

%--------------------------------------------------------------------------
% Step 2: Coil Compression
% -------------------------------------------------------------------------
if iscoil_compression
    fprintf('Coil Compression ... \n');
    nCHA_cc = 14;
    sccmtx = calcSCCMtx(Data_ACS);
    sccmtx_cc = sccmtx(:,1:nCHA_cc);
    compressed_acs_data = CC(Data_ACS,sccmtx_cc);
    compressed_wave_data = CC(rawdata,sccmtx_cc);
else
    compressed_acs_data = Data_ACS;
    compressed_wave_data = rawdata;
    nCHA_cc = size(rawdata,4);
end

%--------------------------------------------------------------------------
% Step 3: ESPIRiT CSM Calibration ~ 6 minutes using 40 cores
% -------------------------------------------------------------------------
fprintf('Bart: ESPIRiT Calibration of full FOV CSM ... \n');
espirit_maps = 1;
compressed_acs_data = zpad(compressed_acs_data,[nRO,Ny_fullFOV,Nz_fullFOV,nCHA_cc]);

bart_command_ecalib = sprintf('ecalib -k%d -m%d -r%d -c0.8 -t0.001 -S',6,espirit_maps,24);
[ecalib_map_fullFOV,eigen_map] = bart(bart_command_ecalib,compressed_acs_data);

weight_mask = eigen_map > 0.8;
weight_mask = repmat(weight_mask, [1,1,1,nCHA_cc,1]);
ecalib_map_fullFOV = ecalib_map_fullFOV .* weight_mask; 


%--------------------------------------------------------------------------
% Step 4: Create multiple-set CSM & PSF for soft-SENSE recon
% -------------------------------------------------------------------------

Nx = size(ecalib_map_fullFOV,1);
Nc = size(ecalib_map_fullFOV,4);

PsfYZ_fullFOV_zpad = zpad(PsfYZ_fullFOV, [Nx_os,Ny_redFOV*2,Nz_redFOV*2]);

PsfYZ_4sets = permute(reshape(circshift(PsfYZ_fullFOV_zpad, [0,Ny_redFOV/2,0]), ...
    [Nx_os,Ny_redFOV,2,Nz_redFOV*2]), [1,2,4,3]);

PsfYZ_4sets = permute(reshape(circshift(PsfYZ_4sets, [0,0,Nz_redFOV/2,Nz_redFOV/2]), ...
    [Nx_os,Ny_redFOV,Nz_redFOV,2,2]), [1,2,3,5,4]);

PsfYZ_4sets = reshape(PsfYZ_4sets,[Nx_os,Ny_redFOV,Nz_redFOV,4]);


CSM_fullFOV_zpad = zpad(ecalib_map_fullFOV, [Nx,Ny_redFOV*2,Nz_redFOV*2,Nc]);

CSM_4sets = permute(reshape(circshift(CSM_fullFOV_zpad, [0,Ny_redFOV/2,0,0]),...
    [Nx,Ny_redFOV,2,Nz_redFOV*2,Nc]), [1,2,4,5,3]);

CSM_4sets = permute(reshape(circshift(CSM_4sets, [0,0,Nz_redFOV/2,0,Nz_redFOV/2]), ...
    [Nx,Ny_redFOV,Nz_redFOV,2,Nc,2]), [1,2,3,5,6,4]);

CSM_4sets = reshape(CSM_4sets,[Nx,Ny_redFOV,Nz_redFOV,Nc,4]);


%--------------------------------------------------------------------------
% Step 5: cg-Soft-SENSE Reconstruction ~ 2 minutes
% -------------------------------------------------------------------------
nIterCG = 30; 

% soft-SENSE reconstruction with 2 sets maps
PsfYZ_4sets = reshape(PsfYZ_4sets,[Nx_os,Ny_redFOV,Nz_redFOV,1,size(PsfYZ_4sets,4)]);
imSize = [Nx,Ny_redFOV,Nz_redFOV,1,size(CSM_4sets,5)];
mask = repmat(sum(abs(compressed_wave_data),4)>0,...
    [1,1,1,size(compressed_wave_data,4)]);
tic
[img,RESVEC] = cgSoftSENSE_wave_v3(compressed_wave_data, ...
     zeros(imSize), CSM_4sets, PsfYZ_4sets, mask, nIterCG);
toc


%--------------------------------------------------------------------------
% Step 6: FOV adjustment to get aliasing-free full FOV image
% -------------------------------------------------------------------------
[Nx,Ny,Nz,Nsets] = size(squeeze(img));

img_fullFOV_Nz = reshape(img,[Nx,Ny,Nz,2,2]);

img_fullFOV_Nz = permute(img_fullFOV_Nz,[1,2,4,3,5]);

img_fullFOV_Nz = flip(img_fullFOV_Nz,5);
tmp = circshift(img_fullFOV_Nz(:,:,:,:), [0,0,0,Nz/2]);
img_fullFOV_NzNy = permute(tmp,[1,4,2,3]);
img_fullFOV_NzNy = flip(img_fullFOV_NzNy,4);
tmp2 = circshift(img_fullFOV_NzNy(:,:,:), [0,0,Ny/2]);
img_final = permute(tmp2,[1,3,2]);

img_fullFOV = crop(img_final,[Nx,Ny_fullFOV,Nz_fullFOV]);
as(img_fullFOV)
