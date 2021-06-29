function [Data_under_redFOV,Data_under_fullFOV,kCalib_fullFOV,kCalib_redFOV,ksp] = ...
CSPI_RedFOV_simulation

fullksp = load('./data/fullksp_3slices.mat');
mask = load('./data/vdpds_mask_new.mat');



ksp = fullksp.ksp3;
ksp = ksp(2:end-1,:,:);
ksp = padarray(ksp,[0,1,0],0);

acssize_y = 24;
acssize_z = 24;

mask_fFOV = mask.vdpds_mask_fFOV;
mask_rFOV = mask.vdpds_mask_rFOV;


Ny_fFOV = size(ksp,1);
Nz_fFOV = size(ksp,2);
Nc = size(ksp,3);

Ny_rFOV = size(mask_rFOV,1);
Nz_rFOV = size(mask_rFOV,2);


% nominal reduced FOV acceleration factors
% i.e. times of fold-over caused by reduced FOV
Rfy = ceil(Ny_fFOV / Ny_rFOV);
Rfz = ceil(Nz_fFOV / Nz_rFOV);


% virtual FOV: multiple of reduced FOV and >= full FOV
Ny_vFOV = Ny_rFOV * Rfy;
Nz_vFOV = Nz_rFOV * Rfz;

Img = ifft2c(ksp);


kCalib_fullFOV = crop(ksp,[acssize_y,acssize_z,Nc]);

Data_under_fullFOV = ksp .* repmat(mask_fFOV,[1,1,Nc]);

% csm: full FOV zpad to virtual FOV
Img_vFOV = zpad(Img,[Ny_vFOV,Nz_vFOV,Nc]);

% shift
shift_amount_y = floor(Ny_vFOV/2) - floor(Ny_rFOV/2);
shift_amount_z = floor(Nz_vFOV/2) - floor(Nz_rFOV/2);
Img_shift = circshift(Img_vFOV,[-shift_amount_y,-shift_amount_z,0]);


% create (Rfy * Rfz) sens maps according to the fold-over process
Img_shift_reshape = reshape(Img_shift,[Ny_rFOV,Rfy,Nz_rFOV,Rfz,Nc]);
Img_shift_reshape = permute(Img_shift_reshape,[1,3,5,2,4]);

tmp = sum(Img_shift_reshape,4);
tmp2 = sum(tmp,5);
ksp_red = fft2c(tmp2);

Data_under_redFOV = ksp_red .* repmat(mask_rFOV,[1,1,Nc]);

kCalib_redFOV = crop(ksp_red,[acssize_y,acssize_y,Nc]);
