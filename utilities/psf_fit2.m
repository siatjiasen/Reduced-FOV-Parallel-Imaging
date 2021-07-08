function [psf_lr,psf_raw,cut_pos_kx,cut_pos_yz] = psf_fit2(psf_raw)
%% - fitting
[Nx,Ny] = size(psf_raw);

lc = round(Ny/4);
Ny_cut = lc*2;
Nx_center = round(Nx/2);
Ny_center = round(Ny/2);

cut_pos_kx = (Nx/2+1-round(Nx/2)):(Nx/2+round(Nx/2));
cut_pos_yz = (Ny_center+1-lc):(Ny_center+lc);

psf_cut = psf_raw(cut_pos_kx,cut_pos_yz);
% [~,Ny_cut] = size(psf_cut);

% unwrap
% psf_cut = unwrap(angle(psf_cut),[],1);
% psf_cut = unwrap(angle(psf_cut),[],2);

lc = floor(size(psf_cut,2)/2);
A = unwrap(angle(psf_cut(:,lc:-1:1)),[],2); 
B = unwrap(angle(psf_cut(:,lc+1:end)),[],2); 
psf_cut = [A(:,lc:-1:1),B];
psf_cut = unwrap(psf_cut,[],1);


% linear fitting
lin = linspace(-1,1,Ny);
lin_cut = lin(cut_pos_yz);

A_mat = cat(2, ones(Ny_cut,1),lin_cut');
A_Mat = cat(2, ones(Ny,1), lin');
coef_y = A_mat \ permute(psf_cut, [2,1]);

% psf-y final estimate
psf_lr = permute(A_Mat * coef_y,[2,1]);
psf_lr = exp(1i * psf_lr);


