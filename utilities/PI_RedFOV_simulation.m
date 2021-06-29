
function [Data_under_redFOV,Data_under_fullFOV,kCalib_fullFOV,kCalib_redFOV,ksp] = ...
PI_RedFOV_simulation(Ry_redFOV,Rz_redFOV, Ry_PI_rFOV,Rz_PI_rFOV,Ry_PI_fFOV,Rz_PI_fFOV)

fullksp = load('./data/fullksp_3slices.mat');

ksp = fullksp.ksp3;
ksp = ksp(2:end-1,:,:);
ksp = padarray(ksp,[0,1,0],0);

acssize_y = 24;
acssize_z = 24;


Ny_fFOV = size(ksp,1);
Nz_fFOV = size(ksp,2);
Nc = size(ksp,3);

Ny_rFOV = 2*round(Ny_fFOV/(2*Ry_redFOV));
Nz_rFOV = 2*round(Nz_fFOV/(2*Rz_redFOV));


% nominal reduced FOV acceleration factors
% i.e. times of fold-over caused by reduced FOV
Rfy = ceil(Ny_fFOV / Ny_rFOV);
Rfz = ceil(Nz_fFOV / Nz_rFOV);


% virtual FOV: multiple of reduced FOV and >= full FOV
Ny_vFOV = Ny_rFOV * Rfy;
Nz_vFOV = Nz_rFOV * Rfz;

Img = ifft2c(ksp);


kCalib_fullFOV = crop(ksp,[acssize_y,acssize_z,Nc]);

Data_under_fullFOV = zeros(size(ksp));

Data_under_fullFOV(1:Ry_PI_fFOV:end,1:Rz_PI_fFOV:end,:) = ...
    ksp(1:Ry_PI_fFOV:end,1:Rz_PI_fFOV:end,:);


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

Data_under_redFOV  = zeros(size(ksp_red));

Data_under_redFOV (1:Ry_PI_rFOV:end,1:Rz_PI_rFOV:end,:) = ...
    ksp_red(1:Ry_PI_rFOV:end,1:Rz_PI_rFOV:end,:);


kCalib_redFOV = crop(ksp_red,[acssize_y,acssize_y,Nc]);

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

