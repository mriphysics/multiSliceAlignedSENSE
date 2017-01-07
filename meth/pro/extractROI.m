function x=extractROI(x,ROI,dim,ext)

%EXTRACTROI3D constructs an image for a given ROI
%   X=EXTRACTROI3D(X,ROI,DIM,EXT) extracts a given ROI from the image or
%   fills the full image
%   X is the input image
%   ROI is the ROI to be extracted / fill the image to
%   DIM is the dimension on which to extract the ROI or fill the full image
%   EXT is a flag that determines whether to extract or fill
%   It returns:
%   X, the output image
%

if ext==1      
    if dim==1
        x=x(ROI(dim,1):ROI(dim,2),:,:,:,:,:,:,:);
    elseif dim==2
        x=x(:,ROI(dim,1):ROI(dim,2),:,:,:,:,:,:);
    elseif dim==3
        x=x(:,:,ROI(dim,1):ROI(dim,2),:,:,:,:,:);
    end
else
    if dim==1
        x=padarray(x,ROI(dim,5),0,'pre');
        x=padarray(x,ROI(dim,6),0,'post');
    elseif dim==2
        x=padarray(x,[0 ROI(dim,5)],0,'pre');
        x=padarray(x,[0 ROI(dim,6)],0,'post');
    elseif dim==3
        x=padarray(x,[0 0 ROI(dim,5)],0,'pre');
        x=padarray(x,[0 0 ROI(dim,6)],0,'post');
    end
end
