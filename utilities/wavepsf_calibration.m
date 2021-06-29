function PsfYZ = wavepsf_calibration(Data_Py,Data_WavePy,Data_Pz,Data_WavePz)

Data_Py = (Data_Py(:,:,:,:,1));
Data_WavePy = (Data_WavePy(:,:,:,:,1));

Data_Pz = (Data_Pz(:,:,:,:,1));
Data_WavePz = (Data_WavePz(:,:,:,:,1));

Data_Py = fftshift(ifft(fftshift(Data_Py,2),[],2),2);
Data_WavePy = fftshift(ifft(fftshift(Data_WavePy,2),[],2),2);
psfy_raw = exp(1i * angle(squeeze(mean(mean(Data_WavePy.*conj(Data_Py),3),4))));
PsfY = psf_fit2(psfy_raw);

Data_Pz = fftshift(ifft(fftshift(Data_Pz,3),[],3),3);
Data_WavePz = fftshift(ifft(fftshift(Data_WavePz,3),[],3),3);
psfz_raw = exp(1i * angle(squeeze(mean(mean(Data_WavePz.*conj(Data_Pz),2),4))));
PsfZ = psf_fit2(psfz_raw);

PsfYZ = repmat(PsfY, [1,1,size(PsfZ,2)]) .* repmat(permute(PsfZ,[1,3,2]),...
    [1,size(PsfY,2),1]); 

end