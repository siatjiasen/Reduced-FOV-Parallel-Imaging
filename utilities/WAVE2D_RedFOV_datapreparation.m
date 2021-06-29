function [DATA_wave_rFOV,DATA_cartesian_rFOV,DATA_wave_fFOV,DATA_cartesian_fFOV,kCalib_fullFOV,kCalib_redFOV,PsfY_fFOV,PsfY_rFOV] = ...
WAVE2D_RedFOV_datapreparation(Ry_rFOV,Ry_fFOV)


load('./data/Figure4_2D.mat');

Nx_os = size(DATA_wave_fFOV,1);
Nx = size(DATA_cartesian_fFOV,1);

Ny_fFOV = size(DATA_wave_fFOV,2);
Ny_rFOV = size(DATA_wave_rFOV,2);

Nc = size(DATA_wave_rFOV,3);


% =========================================================================
% Generate full FOV & reduced FOV calibration data
% =========================================================================

kCalib_redFOV = crop(DATA_cartesian_rFOV, [24,24,Nc]);  % reduced FOV calibration data
kCalib_fullFOV = crop(DATA_cartesian_fFOV, [24,24,Nc]); %  full FOV calibration data

% =========================================================================
% Generate under-sampled reduced FOV wave data
% =========================================================================
mask_wave = zeros(Nx_os,Ny_rFOV);
mask_wave(:,end/2+1:-Ry_rFOV:1) = 1;
mask_wave(:,end/2+1:Ry_rFOV:end) = 1;
samp_mat_wave = repmat(mask_wave,[1,1,Nc]);
DATA_wave_rFOV = DATA_wave_rFOV .* samp_mat_wave;


mask = zeros(Nx,Ny_rFOV);
mask(:,end/2+1:-Ry_rFOV:1) = 1;
mask(:,end/2+1:Ry_rFOV:end) = 1;
samp_mat = repmat(mask,[1,1,Nc]);
DATA_cartesian_rFOV = DATA_cartesian_rFOV .* samp_mat;

% =========================================================================
% Generate under-sampled full FOV wave data
% =========================================================================

mask_wave = zeros(Nx_os,Ny_fFOV);
mask_wave(:,end/2+1:-Ry_fFOV:1) = 1;
mask_wave(:,end/2+1:Ry_fFOV:end) = 1;
samp_mat_wave = repmat(mask_wave,[1,1,Nc]);
DATA_wave_fFOV = DATA_wave_fFOV .* samp_mat_wave;


mask = zeros(Nx,Ny_fFOV);
mask(:,end/2+1:-Ry_fFOV:1) = 1;
mask(:,end/2+1:Ry_fFOV:end) = 1;
samp_mat = repmat(mask,[1,1,Nc]);
DATA_cartesian_fFOV = DATA_cartesian_fFOV .* samp_mat;



% Im_fFOV = ifft2c(DATA_fFOV);
% Im_fFOV_zpadx = zpad(Im_fFOV,[size(PsfY_fFOV,1),size(Im_fFOV,2),size(Im_fFOV,3)]);
% DATA_wave_fFOV = fftc(Im_fFOV_zpadx,1) .* repmat(PsfY_fFOV,[1,1,size(Im_fFOV,3)]);
% DATA_wave_fFOV = fftc(DATA_wave_fFOV,2);