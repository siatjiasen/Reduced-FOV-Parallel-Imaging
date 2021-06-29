function [res,RESVEC,LSVEC,FLAG] = cgSoftSENSE_wave_v3(kData, x0, csm, PsfYZ ,mask, nIterCG, alpha)
% conjugate gradient (cg) algo iteration to solve wave recon
%          || K * F_yz * Psf * F_x * C * v - m ||_2 + lambda * R (v)
%     where, m - wave encoded samples;    v - image;  C - 3D CSM
%     K - discrete sampling; F_yz - 2D Fourier Tranform along y and z
%     Psf - 3D PSF in hybrid domain (kx,y,z);  F_x - 1D FT along x
%    
% Forward
% v (im) -> C*im  -> zpad(C*im) along RO -> fft1(kx) -> PsfYZ -> fft2(Kyz)
% Adjoint
% ifft2(k_yz) -> conj(PsfYZ) -> ifft1(RO) -> crop() along RO -> conj(C) -> k
%
% Authored by Jia Sen 

nIterSplit = 1;
nInnerCG = nIterCG;

if nargin < 7
   alpha = 0.01;
end

imSize = size(x0);
if length(imSize) == 2
    imSize = [imSize,1,1,1];
end
if length(imSize) == 3
    imSize = [imSize,1,1];
end
dataSize = size(kData);

res = x0(:);

for n=1:nIterSplit
    b = [kData(:); sqrt(alpha)*res(:)];
    
    [res,FLAG,~,~,RESVEC,LSVEC] = lsqr(@(x,tflag)afun(x,csm,PsfYZ,mask,dataSize,imSize,alpha,tflag),...
        b, [], nInnerCG,[],[], res(:));
    res = reshape(res,imSize);
    
end
 

function [y, tflag] = afun(x, csm, PsfYZ, mask, dataSize, imSize,  alpha, tflag)
    
if strcmp(tflag,'transp')
    % Adjoint
    % ifft2(k_yz) -> conj(PsfYZ) -> ifft1(RO) -> crop() along RO -> conj(S) -> k
    
    % Input: multiple coil kspace
    y = reshape(x(1:prod(dataSize)),dataSize);
    xtmp = x(prod(dataSize)+1:end);
    
    % Adjoint of Subsampling operator
    y_ = y.*mask;

    % Adjoint PE 1& 2 FFT operator
    im = ifftc(ifftc(y_,2),3);
    
    % Adjoint Soft-SENSE Wave PSF (Multiple-Set) Operator
    im_wave = bsxfun(@times,conj(PsfYZ),im);

    % Adjoint RO FFT operator
    im = ifftc(im_wave,1);
    
    % Adjoint RO Resize operator
    if imSize(5) == 1 || size(PsfYZ,5) == 1
        x = crop(im,[imSize(1),dataSize(2),dataSize(3),dataSize(4)]);
    else
        x = crop(im,[imSize(1),dataSize(2),dataSize(3),dataSize(4),imSize(5)]);
    end
   
    % Adjoint Soft-SENSE Sensitivity MAP (Multiple-Set) operator
    x = sum(bsxfun(@times,conj(csm),x),4);
    
    % Output - multiple-set images
    y = x(:) + sqrt(alpha) * xtmp(:);
    
    
else
    % Forward
    % im -> S*im  -> zpad(S*im) along RO -> fft1(RO) -> PsfYZ -> fft2(Kyz)
    
    % Input - multiple-set images (the last dim corresponds to ESPIRiT sets)
    x = reshape(x,[imSize(1),imSize(2),imSize(3),1,imSize(5)]);
    
    % Soft-SENSE Sensitivity MAP (Multiple-Set) operator
    x_ = bsxfun(@times,x,csm);
    
    % RO Resize operator
    x_zpadRO = zpad(x_,[dataSize(1),imSize(2),imSize(3),dataSize(4),imSize(5)]);
    
    % RO FFT operator
    k_ = fftc(x_zpadRO,1);

    % Soft-SENSE Wave PSF (Multiple-Set) operator
    PsfYZ = reshape(PsfYZ,[dataSize(1),imSize(2),imSize(3),1,size(PsfYZ,5)]);
    k_PsfYZ = sum(bsxfun(@times, k_, PsfYZ),5);
    
    % PE1&2 FFT operator
    y = fftc(fftc(k_PsfYZ,2),3);
    
    % Subsampling operator
    y = y .* mask;
    
    % Output - multiple-coil kspace
    y = [y(:); sqrt(alpha) * x(:)];
    
end