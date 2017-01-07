function ROI=computeROI(W,ext,strict)

%COMPUTEROI computes the ROI for reconstruction
%   ROI=COMPUTEROI(W,EXT,STRICT) computes the spatial ROI for
%   reconstruction corresponding to a given spatial mask
%   W is the mask of the structure of interest
%   EXT is a range of security of the ROI with respect to the mask
%   STRICT is a flag to select whether the resulting ROI along a given
%   direction has to be kept centered with respect to the acquired FOV
%   It returns:
%   ROI, ranges of the ROI (columns 1 and 2), size of the original image
%   (column 3) and other information as specified below (remaining columns)
%

ND=ndims(W);
if ~exist('strict','var');strict=zeros(1,ND);end
N=size(W);
ROI=zeros(ND,6);
ROI(:,1)=ones(ND,1);%Inferior limit of the ROI
ROI(:,2)=N';%Superior limit of the ROI
ROI(:,3)=N';%Size of the original image
a(1)=1;
a(2)=-1;
for l=1:2
    for nd=1:ND
        Wu=shiftdim(W,nd-1);
        while isempty(find(Wu(ROI(nd,l),:,:)~=0,1))
            ROI(nd,l)=ROI(nd,l)+a(l);            
        end
    end
end
           
for nd=1:ND
    ROI(nd,1)=max([ROI(nd,1)-ext 1]);     
    ROI(nd,2)=min([ROI(nd,2)+ext ROI(nd,3)]);  
    ROI(nd,4)=ROI(nd,2)-ROI(nd,1)+1;%Number of elements in the ROI
    ROI(nd,5)=ROI(nd,1)-1;%Number of discarded elements at the left
    ROI(nd,6)=ROI(nd,3)-ROI(nd,2);%Number of discarded elements at the right
    if strict(nd)
        [~,indSm]=min(ROI(nd,5:6));
        if indSm==1
            ROI(nd,4)=ROI(nd,4)+ROI(nd,6)-ROI(nd,5);
            ROI(nd,2)=ROI(nd,2)+ROI(nd,6)-ROI(nd,5);
            ROI(nd,6)=ROI(nd,5);
        else
            ROI(nd,4)=ROI(nd,4)+ROI(nd,5)-ROI(nd,6);
            ROI(nd,1)=ROI(nd,1)+ROI(nd,5)-ROI(nd,6);
            ROI(nd,5)=ROI(nd,6);
        end
    end            
end
