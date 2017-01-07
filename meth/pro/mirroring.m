function x=mirroring(x,mirror,di,ty)

%MIRRORING mirrors a given image so as to guarantee periodicity
%   X=MIRRORING(X,MIRROR,DI,TY) mirrors an image along certain directions
%   Inputs:
%   X is the image to be mirrored
%   MIRROR is a flag to select whether a given direction has to be mirrored
%   DI is the direction of mirroring: 1->mirror / 0->demirror
%   TY is the type of mirroring: 0->unsymmetric / 1->symmetric (defaults to
%   0)
%   It returns:
%   X, the mirrored image
%

if ~exist('ty','var');ty=0;end

N=size(x);
nDimsOu=length(mirror);
N(end+1:nDimsOu)=1;
nDimsIn=length(N);
if nDimsOu>nDimsIn
    error('Mirroring dimensionality is larger than image dimensionality');
end
if nDimsIn>12
    error('Mirroring only works for arrays of up to 12 dimensions');
end

for m=1:nDimsOu       
    if mirror(m)       
        perm=1:nDimsIn;perm(1)=m;perm(m)=1;                
        x=permute(x,perm);
        if di            
            if ~ty
                x=[x;x(end:-1:1,:,:,:,:,:,:,:,:,:,:,:)];
            else
                x=[x(end:-1:1,:,:,:,:,:,:,:,:,:,:,:);x;x(end:-1:1,:,:,:,:,:,:,:,:,:,:,:)];
            end
        else
            if ~ty
                x=x(1:end/2,:,:,:,:,:,:,:,:,:,:,:);
            else
                x=x(end/3+1:2*end/3,:,:,:,:,:,:,:,:,:,:,:);
            end
        end
        x=permute(x,perm);
    end
end
