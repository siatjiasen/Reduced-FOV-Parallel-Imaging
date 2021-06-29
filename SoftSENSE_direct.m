function [Img_recon,G_map] = SoftSENSE_direct(Data_under,CSM,im,Ry,Rz)
% Aliasing-free Soft-SENSE reconstruction of 
%   uniformly downsampled parallel imaging with full or reduced FOV imaging
%   using single or multiple-set CSM estimated from full FOV calibration scan
% Input 
%    Data_under - acquired multi-coil k-space (full or reduced FOV kspace)
%    CSM - single or multiple-set coil sensitivity maps
%    im  - initial value of reconstructed image, zeros(im_size)
%    Ry - parallel imaging acceleration factor along Phase Encoding 1
%    (Integer)
%    Rz - parallel imaging acceleration factor along Phase Encoding 2
%    (Integer)
% Output
%    Img_recon - Full FOV image
%    G-factor  - Full FOV g-factor
%
% Authored by Jia Sen

% acquired kspace size
Ny = size(Data_under,1);
Nz = size(Data_under,2);
Nc = size(Data_under,3);

% image matrix size
Ny_fFOV = size(im,1);
Nz_fFOV = size(im,2);

% nominal reduced FOV acceleration factors
% i.e. times of fold-over caused by reduced FOV
Rfy = ceil(Ny_fFOV / Ny);
Rfz = ceil(Nz_fFOV / Nz);

% virtual FOV: multiple of reduced FOV and >= full FOV
Ny_vFOV = Ny * Rfy;
Nz_vFOV = Nz * Rfz;

% shift
shift_amount_y = floor(Ny_vFOV/2) - floor(Ny/2);
shift_amount_z = floor(Nz_vFOV/2) - floor(Nz/2);

% under-sampled kspace transformed to aliased image
Img_alias = ifft2c(Data_under);

% cut the first (Ny/Ry,Nz/Rz) region
Img_alias = Img_alias(1:Ny/Ry,1:Nz/Rz,:);

% unfold aliased voxels by matrix inverse
Img_recon = zeros(Ny,Nz,Rfy,Rfz);
G_map = zeros(Ny,Nz,Rfy,Rfz);

tic;
for iz = 1:Nz/Rz
    for iy = 1:Ny/Ry
        idx_y = iy:Ny/Ry:Ny; 
        idx_z = iz:Nz/Rz:Nz; 
        b = squeeze(Img_alias(iy,iz,:));
        A = permute(CSM(idx_y,idx_z,:,:,:),[1,2,4,5,3]);
        A = permute(reshape(A,[Ry*Rz*Rfy*Rfz,Nc]),[2,1]);
        x = pinv(A) * b; % matrix pseudo inverse
        Img_recon(idx_y,idx_z,:,:) = reshape(x,[Ry,Rz,Rfy,Rfz]);
        
        % calculate g-factors
        g_rho = sqrt(real(diag(pinv(A'*A)).*diag(A'*A)));
        G_map(idx_y,idx_z,:,:) = reshape(g_rho,[Ry,Rz,Rfy,Rfz]);
    end
end
toc;

% joint image: virtual FOV
Img_recon = reshape(permute(Img_recon,[1,3,2,4]),[Ny*Rfy,Nz*Rfz]);
G_map = reshape(permute(G_map,[1,3,2,4]),[Ny*Rfy,Nz*Rfz]);

% shift back
Img_recon = circshift(Img_recon,[shift_amount_y,shift_amount_z]);
G_map = circshift(G_map,[shift_amount_y,shift_amount_z]);

% virtual FOV cut to full FOV
Img_recon = crop(Img_recon,[Ny_fFOV,Nz_fFOV]);
G_map = crop(G_map,[Ny_fFOV,Nz_fFOV]);

end