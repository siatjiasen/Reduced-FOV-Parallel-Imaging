%-------------------------------------------------------------------------%
% Code for reproducing *** Figure 4 *** for
%     "Aliasing-free Reduced Field-of-View Parallel Imaging"
% Note: Run "WAVE2D_RedFOV_datapreparation" first for data preparation
%-------------------------------------------------------------------------%
clc;clear;close all

addpath(strcat(pwd,'\data'))
addpath(genpath(strcat(pwd,'\utilities')))

do_waverecon_using_softsense = 1;  % 0-SPIRiT model, which needs 2-3x more computational time

% =========================================================================
% Load Data
% =========================================================================
Ry_rFOV = 3;
Ry_fFOV = 5;

[DATA_wave_redFOV,DATA_cartesian_redFOV,DATA_wave_fullFOV,DATA_cartesian_fullFOV,...
    kCalib_fullFOV,kCalib_redFOV,PsfY_fullFOV,PsfY_redFOV] = ...
    WAVE2D_RedFOV_datapreparation(Ry_rFOV,Ry_fFOV);

% Get data size
[Nx,Ny_fullFOV,Nc] = size(DATA_cartesian_fullFOV);
Ny_redFOV = size(DATA_wave_redFOV,2);        % reduced FOV phase encoding size
Nx_os = size(DATA_wave_redFOV,1);       % oversampled kx size due to wave encoding


if do_waverecon_using_softsense
    % =========================================================================
    % Estimate full FOV (1 set) and reduced FOV (2 sets) CSM
    % =========================================================================
    disp('ESPIRiT calibration of full FOV CSM and multiple-set reduced FOV CSM')
    kSize = [6,6]; eigThresh_1 = 0.02; eigThresh_2 = 0.95;
    imSize = [Nx,Ny_fullFOV]; Nmaps = 1;
    ecalib_map_fullFOV = espirit_calib2d(kCalib_fullFOV,imSize,kSize,eigThresh_1,eigThresh_2,Nmaps);
    
    kSize = [6,6]; eigThresh_1 = 0.02; eigThresh_2 = 0.9;
    imSize = [Nx,Ny_redFOV]; Nmaps = 2;
    ecalib_map_redFOV = espirit_calib2d(kCalib_redFOV,imSize,kSize,eigThresh_1,eigThresh_2,Nmaps);
    ecalib_map_redFOV = flip(ecalib_map_redFOV,4);
    
    % =========================================================================
    % Generate 2-set CSM and PSF from 1 set of full FOV maps and PSF
    % =========================================================================
    disp('Generate multiple-set full FOV CSM and PSF for wave-soft-SENSE')
    ecalib_maps_fullFOV_zpad = zpad(ecalib_map_fullFOV, [Nx,2*Ny_redFOV,Nc]);
    ecalib_map_fullFOV_2sets = flip(reshape(circshift(ecalib_maps_fullFOV_zpad, ...
        [0,Ny_redFOV/2,0]), [Nx,Ny_redFOV,2,Nc]),3);
    
    PsfY_fullFOV_zpad = zpad(PsfY_fullFOV, [Nx_os,2*Ny_redFOV]);
    PsfY_fullFOV_2sets = reshape(circshift(PsfY_fullFOV_zpad, ...
        [0,Ny_redFOV/2,0]),[Nx_os,Ny_redFOV,2]);
    PsfY_fullFOV_2sets = flip(PsfY_fullFOV_2sets,3);
    PsfY_fullFOV_2sets(abs(PsfY_fullFOV_2sets) == 0) = 1;    % important for convergence rate
    
    
    % =========================================================================
    % Reconstruction: using different choices of CSM and PSF
    % =========================================================================
    disp('Conjugate Gradient Solver')
    nIterCG = 60;
    lambda = 0.003;
    imSize = [Nx,Ny_redFOV,1,1,2];
    compressed_wave_data = reshape(DATA_wave_redFOV,[Nx_os,Ny_redFOV,1,Nc]);
    mask = repmat(sum(abs(compressed_wave_data),4)>0,...
        [1,1,1,size(compressed_wave_data,4)]);
    
    % 2-set full FOV CSM and 2-set full FOV wave-PSF
    disp('1. Aliasing-free recon using 2-set full FOV CSM and 2-set full FOV wave-PSF')
    CSM_2sets = reshape(permute(ecalib_map_fullFOV_2sets,[1,2,4,3]),[Nx,Ny_redFOV,1,Nc,2]);
    PsfYZ_2sets = reshape(PsfY_fullFOV_2sets,[Nx_os,Ny_redFOV,1,1,2]);
    
    tic
    img_CSM2_PSF2 = cgSoftSENSE_wave_v3(compressed_wave_data, ...
        zeros(imSize), CSM_2sets, PsfYZ_2sets, mask, nIterCG, lambda);
    toc
    
    as(squeeze(img_CSM2_PSF2))
    
    
    % 2-set reduced FOV CSM and 1-set reduced FOV wave-PSF
    disp('2. using 2-set ESPIRiT CSM and 1-set reduced FOV wave-PSF leads to severe aliasing')
    CSM_2sets = reshape(ecalib_map_redFOV,[Nx,Ny_redFOV,1,Nc,2]);
    PsfYZ = reshape(PsfY_redFOV,[Nx_os,Ny_redFOV,1,1,1]);
    
    tic
    img_CSMredFOV_PSFredFOV = cgSoftSENSE_wave_v3(compressed_wave_data, ...
        zeros(imSize), CSM_2sets, PsfYZ, mask, nIterCG, lambda);
    toc
    
    as(squeeze(img_CSMredFOV_PSFredFOV))
    
    % 2-set reduced FOV CSM and 2-set full FOV wave-PSF
    disp('3. using 2-set ESPIRiT CSM and 2-set full FOV wave-PSF leads to severe aliasing')
    CSM_2sets = reshape(ecalib_map_redFOV,[Nx,Ny_redFOV,1,Nc,2]);
    PsfYZ_2sets = reshape(PsfY_fullFOV_2sets,[Nx_os,Ny_redFOV,1,1,2]);
    
    tic
    img_CSMredFOV_PSFfullFOV = cgSoftSENSE_wave_v3(compressed_wave_data, ...
        zeros(imSize), CSM_2sets, PsfYZ_2sets, mask, nIterCG,lambda);
    toc
    
    as(squeeze(img_CSMredFOV_PSFfullFOV))
    
    % 2-set full FOV CSM and 1-set reduced FOV wave-PSF
    disp('4. using 2-set full FOV CSM and 1-set reduced FOV wave-PSF leads to severe aliasing')
    CSM_2sets = reshape(permute(ecalib_map_fullFOV_2sets,[1,2,4,3]),[Nx,Ny_redFOV,1,Nc,2]);
    PsfYZ = reshape(PsfY_redFOV,[Nx_os,Ny_redFOV,1,1,1]);
    
    tic
    img_CSMfullFOV_PSFredFOV = cgSoftSENSE_wave_v3(compressed_wave_data, ...
        zeros(imSize), CSM_2sets, PsfYZ, mask, nIterCG, lambda);
    toc
    
    as(squeeze(img_CSMfullFOV_PSFredFOV))
    
    % Cartesian reduced FOV reconstruction using 2-set full FOV CSM
    disp('5. Reduced FOV Cartesian Imaging using 2-set full FOV CSM - Noisy withou Wave encoding')
    compressed_cartesian_data = reshape(DATA_cartesian_redFOV,[Nx,Ny_redFOV,1,Nc]);
    mask = repmat(sum(abs(compressed_cartesian_data),4)>0,...
        [1,1,1,size(compressed_cartesian_data,4)]);
    
    CSM_2sets = reshape(permute(ecalib_map_fullFOV_2sets,[1,2,4,3]),[Nx,Ny_redFOV,1,Nc,2]);

    tic
    img_CSM2_cartesian = cgSoftSENSE_cartesian_v3(compressed_cartesian_data, ...
        zeros(imSize), CSM_2sets, mask, nIterCG, 0.002);
    toc
    as(squeeze(img_CSM2_cartesian))
    
else
    % =========================================================================
    % SPIRiT Reconstruction for Wave-Encoded data
    % =========================================================================
    disp('performing calibration for SPIRiT');
    kSize = [5,5];      % SPIRiT kernel size
    
    tic;
    kernel_fullFOV = calibSPIRiT(kCalib_fullFOV, kSize, Nc, 0.05);
    GOP_fullFOV = SPIRiT(kernel_fullFOV, 'image',[Nx,Ny_fullFOV]);
    toc;
    
    disp('Conjugate Gradient Solver for Wave-SPIRiT');
    nIterCG = 60;
    PsfY_1 = reshape(PsfY_fullFOV,[Nx_os,Ny_fullFOV,1,1]);
    x0 = zeros(Nx,Ny_fullFOV,1,Nc);
    y = reshape(DATA_wave_redFOV,[Nx_os,Ny_redFOV,1,Nc]);
    lambda = 1;
    tic
    [res_red2full_wavespirit, RESVEC] = cgSPIRiT_wave_img_v6(y, GOP_fullFOV, PsfY_1, nIterCG, lambda, x0);
    toc
    im_red2full_wavespirit = sos(res_red2full_wavespirit);
    as(im_red2full_wavespirit)

    y = reshape(DATA_cartesian_redFOV,[Nx,Ny_redFOV,1,Nc]);
    lambda = 1;
    tic
    [res_red2full_spirit, RESVEC] = cgSPIRiT_cartesian_img_v6(y, GOP_fullFOV, nIterCG, lambda, x0);
    toc
    im_red2full_spirit = sos(res_red2full_spirit);
    as(im_red2full_spirit)
    
end


