function [res,RESVEC,LSVEC,FLAG] = cgSoftSENSE_cartesian_v3(kData, x0, csm,mask, nIterCG, alpha)
% Soft-SENSE recon by conjugate gradient (cg) algo to solving:
%          || K * F_xyz * C * v - m ||_2 + lambda * R (v)
%     where, m - wave encoded samples;    v - image;  C - 3D CSM
%     K - discrete sampling; F_xyz - 3D FT along x
%    
% Forward
% v (im) -> C*im  ->  fft3(Kxyz)
% Adjoint
% ifft3(kxyz) -> conj(C) -> k (kspace)
%
% Authored by Jia Sen 


nIterSplit = 1;        
nInnerCG = nIterCG;

if nargin < 6
    alpha = 0.001;
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

    [res,FLAG,~,~,RESVEC,LSVEC] = lsqr(@(x,tflag)afun(x,csm,mask,dataSize,imSize,alpha,tflag),...
        b, [], nInnerCG,[],[], res(:));
    res = reshape(res,imSize);   
    
end
 

function [y, tflag] = afun(x, csm, mask, dataSize, imSize, alpha, tflag)
    
if strcmp(tflag,'transp')
    % Adjoint : ifft3(k_xyz) -> conj(S) -> k
    
    % Input: multiple coil kspace  (reduced FOV)
    y = reshape(x(1:prod(dataSize)),dataSize);
    xtmp = x(prod(dataSize)+1:end);
    
    % Adjoint of Subsampling operator
    y_ = y.*mask;

    % Adjoint 3D FFT operator
    im = ifftc(ifftc(ifftc(y_,3),2),1);
 
    % Adjoint Soft-SENSE Sensitivity MAP operator:  unaliasing to multiple-set images
    x = sum(bsxfun(@times,conj(csm),im),4); 
  
    % Output - multiple-set images (the last dim corresponds to softSENSE sets)
    y = x(:) + sqrt(alpha) * xtmp(:);
    
else
    % Forward  :  im -> S*im  -> fft3(k_xyz)
    
    % Input - multiple-set images (the last dim corresponds to softSENSE sets)
    x = reshape(x,[imSize(1),imSize(2),imSize(3),1,imSize(5)]);
    
    % Soft-SENSE Sensitivity MAP operator: sum multiple-set images to aliased image
    x_ = sum(bsxfun(@times,csm,x),5);
   
    % 3D FFT operator
    k_ = fftc(fftc(fftc(x_,1),2),3);
    
    % Subsampling operator
    y = k_.* mask;
    
    % Output - multiple-coil kspace (reduced FOV)
    y = [y(:); sqrt(alpha) * x(:)];
    
end