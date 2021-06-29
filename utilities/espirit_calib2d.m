function ecalib_map = espirit_calib2d(kCalib,imSize,kSize,eigThresh_1,eigThresh_2,Nmaps)
tic;

Ny = imSize(1);
Nz = imSize(2);
Nc = size(kCalib,3);

[K,S] = dat2Kernel(kCalib, kSize);
idx = find(S >= S(1)*eigThresh_1, 1, 'last');
[M,W] = kernelEig(K(:,:,:,1:idx), [Ny,Nz]);

weight_mask = W(:,:,end-Nmaps+1:end) > eigThresh_2;
weight_mask = repmat(permute(weight_mask,[1,2,4,3]), [1,1,Nc,1]);
ecalib_map = M(:,:,:,end-Nmaps+1:end) .* weight_mask; 

toc;

end