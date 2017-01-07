function writeNIFTI(x,Phi,pathO,ima,suffix)

%WRITENIFTI writes a nifti file with the reconstructed volume
%   WRITENIFTI(X,PHI,PATHO,IMA,SUFFIX) saves reconstructions in .nii
%   format
%   X is the image to be saved
%   PHI is the image spacing and orientation
%   PATH is the folder where to save the data
%   IMA is the file where to save the data
%   SUFFIX is the suffix to be added to the filename (defaults to '')
%

if ~exist('suffix','var');suffix='';end

MT=Phi;
MS=sqrt(sum(MT(1:3,1:3).^2,1));%Might present problems in cases of oblique acquisitions...           
niftiIm=make_nii(single(x),MS);
MTT=[-1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1];    
MT=MTT*MT;
MT(1:3,4)=MT(1:3,4)+MT(1:3,3)/MS(3)+MT(1:3,2)/MS(2)+MT(1:3,1)/MS(1);
niftiIm.hdr.hist.srow_x=MT(1,:);
niftiIm.hdr.hist.srow_y=MT(2,:);
niftiIm.hdr.hist.srow_z=MT(3,:);
niftiIm.hdr.hist.sform_code=1;
save_nii(niftiIm,sprintf('%s/%s%s.nii',pathO,ima,suffix));
