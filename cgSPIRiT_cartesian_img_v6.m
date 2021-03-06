function [res, LSVEC] = cgSPIRiT_cartesian_img_v6(y, GOP, nIterCG, lambda, x0)
% Cartesian SPIRiT recon solved by conjugate gradient (cg) algorithm
%      lambda * ||(G - I) * x ||_2  + || D * F * A * x - y ||_2
%      SPIRiT Calibration Consistency + Data Consistency
% where,
%     x - unknown multi-coil images (full FOV), x0 - initial value
%     y - acquired multi-coil k-space (full or reduced FOV kspace)
%     G - SPIRiT kernel in image space, GOP
%     F - 3D Fourier Tranform along x, y and z
%     A - Image Aliasing operator due to reduced FOV imaging
% Supporting 2D (kx,ky,1) or (1,ky,kz) and 3D (kx,ky,kz) 
%     Reduced FOV kspace acquisition, Reduced FOV image reconstruction
%     Full    FOV kspace acquisition, Full    FOV image reconstruction
%  ** Reduced FOV kspace acquisition, Full    FOV image reconstruction **
%
% Authored by Jia Sen

if nargin < 4
    lambda = 0.1;
end

if nargin < 5
    x0 = GOP*ifft2c(y);
end


[Nx,Ny_fullFOV,Nz_fullFOV,Nc] = size(x0); % Image size to be estimated (full FOV)
[Nx_os,Ny_redFOV,Nz_redFOV,~] = size(y);     % kspace size acquired

idx_acq = abs(y) > 0;
yy = [y(idx_acq); zeros(size(x0(:)))];

res = x0(:);
for iter_outer = 1:1:1
    
    [tmpres,~,~,~,~,LSVEC] = lsqr(@aprod,yy,1e-6,nIterCG,[],[],res(:),...
        GOP,Nx,Nx_os,Ny_fullFOV,Nz_fullFOV,Nc,Ny_redFOV,Nz_redFOV,idx_acq,lambda);
    res = reshape(tmpres,[Nx,Ny_fullFOV,Nz_fullFOV,Nc]);
    
end

end


function [res,tflag] = aprod(x,GOP,Nx,Nx_os,Ny_fullFOV,Nz_fullFOV,Nc,Ny_redFOV,Nz_redFOV,idx_acq,lambda,tflag)

if strcmp(tflag,'transp')
    % Input: acquired reduced FOV k-space data
    
    %%% Adjoint of Acquired Data Consistency  Fxyz^H * D^H * y = x
    % D^H * y
    k_rFOV = zeros(Nx_os,Ny_redFOV,Nz_redFOV,Nc);
    k_rFOV(idx_acq) = x(1:end-Nx*Ny_fullFOV*Nz_fullFOV*Nc);
    
    % Fyz^H * D^H * y
    im_rFOV = ifftc(ifftc(ifftc(k_rFOV,1),2),3);
    
    % x = A^H * F^H * D^H * y
    if Ny_redFOV < Ny_fullFOV || Nz_redFOV < Nz_fullFOV
        tmpy = general_aliasing(im_rFOV,Ny_fullFOV,Nz_fullFOV,Nc,Ny_redFOV,Nz_redFOV,'transp');
    else
        tmpy = im_rFOV;
    end
    
    %%% Adjoint of SPIRiT Calibration Consistency (G-I)^H * x = 0
    im_fFOV = reshape(x(end-Nx*Ny_fullFOV*Nz_fullFOV*Nc+1:end),[Nx,Ny_fullFOV,Nz_fullFOV,Nc]);
    tmpx = GOP' * squeeze(im_fFOV);
    
    % Output: full FOV multi-coil images
    res = tmpy(:) + lambda * tmpx(:);
    
else
    
    % Input: full FOV multi-coil images
    im_fFOV = reshape(x,[Nx,Ny_fullFOV,Nz_fullFOV,Nc]);
    
    %%% Acquired Data Consistency: D * Fxyz * x = y
  
    % From full FOV images to aliased reduced FOV images: A * x
    if Ny_redFOV < Ny_fullFOV || Nz_redFOV < Nz_fullFOV
        im_rFOV = general_aliasing(im_fFOV,Ny_fullFOV,Nz_fullFOV,Nc,Ny_redFOV,Nz_redFOV,'notransp');
    else
        im_rFOV = im_fFOV;
    end
    
    % PE1 & PE2 FFT operator
    ksp = fftc(fftc(fftc(im_rFOV,1),2),3);
    
    % Acquired Data Consistency: D * F * A * x
    y_acq = ksp(idx_acq);
    
    %%% SPIRiT Calibration Consistency ("image"): (G-I)*x = 0
    im_fFOV_G = GOP * squeeze(im_fFOV);
    
    % Output: acquired reduced FOV kspace
    res = [y_acq(:); lambda * im_fFOV_G(:)];
    
end
end


function outp = general_aliasing(inp,Ny_fullFOV,Nz_fullFOV,Nc,Ny_redFOV,Nz_redFOV,tflag)

if strcmp(tflag,'transp')
    % Input: aliased reduced FOV multi-coil images
    % Output: full FOV multi-coil images
    
    if Ny_redFOV < Ny_fullFOV
        outp = repmat(inp, [1,3,1,1]);
        outp = crop(outp, [size(inp,1),Ny_fullFOV,size(inp,3),Nc]);
    else
        outp = inp;
    end
    if Nz_redFOV < Nz_fullFOV
        outp = repmat(outp, [1,1,3,1]);
        outp = crop(outp, [size(inp,1),Ny_fullFOV,Nz_fullFOV,Nc]);
    end
    
else
    % Input:  full FOV multi-coil images
    % Output: aliased reduced FOV multi-coil images
    
    if Ny_redFOV < Ny_fullFOV
        inp_zpad = zpad(inp,[size(inp,1),3*Ny_redFOV,size(inp,3),Nc]);
        outp = inp_zpad(:,1:Ny_redFOV,:,:) + inp_zpad(:,Ny_redFOV+1:2*Ny_redFOV,:,:) + ...
            inp_zpad(:,2*Ny_redFOV+1:3*Ny_redFOV,:,:);
    else
        outp = inp;
    end
    if Nz_redFOV < Nz_fullFOV
        outp_zpad = zpad(outp,[size(outp,1),size(outp,2),3*Nz_redFOV,Nc]);
        outp = outp_zpad(:,:,1:Nz_redFOV,:) + outp_zpad(:,:,Nz_redFOV+1:2*Nz_redFOV,:) + ...
            outp_zpad(:,:,2*Nz_redFOV+1:3*Nz_redFOV,:);
    end
    
end
end
