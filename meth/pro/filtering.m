function x=filtering(x,H,gpu)

% FILTERING filters a given image
%   X=FILTERING(X,H,GPU) filters an image
%   X is the image to be filtered
%   H is the filter 
%   GPU is a flag that determines whether to use gpu processing (defaults
%   to 0)
%   It returns:
%   X, the filtered image
%

if ~exist('gpu','var');gpu=0;end

NH=size(H);
nDimsH=ndims(H);

for m=1:nDimsH
    if NH(m)~=1
        x=fftGPU(x,m,gpu);
    end
end
x=bsxfun(@times,x,H);
for m=1:nDimsH
    if NH(m)~=1
        x=ifftGPU(x,m,gpu);
    end
end
