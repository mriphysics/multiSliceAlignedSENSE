function x=gibbsRingingFilter(x,NDims,gibbsRing,gpu)

% GIBBSRINGINGFILTER applies Gibbs ringing filter
%   X=GIBBSRINGINGFILTER(X,NDIMS,GIBBSRING,GPU) applies a Gibbs filter
%   based on a Tukey window
%   X is the image to be filtered
%   NDIMS are the number of dimensions to filter
%   GIBBSRING is the  Gibbs ringing correction factor. Scalar: 0 -> no 
%   correction / 1 -> Hann window. Vectorial: anisotropic
%   GPU is a flag that determines whether to use gpu processing (defaults
%   to 0)
%   It returns:
%   X, the filtered image
%

if ~exist('gpu','var');gpu=0;end

N=size(x);
NGR=length(gibbsRing);
if NGR==1    
    NNDims=N(1:NDims);
    kk=single(zeros(prod(NNDims),NDims));
    if gpu>0
        kk=gpuArray(kk);
    end
    kGrid=generateGrid(NNDims,1,NNDims,[],1,gpu);
    for m=1:NDims    
        rep=NNDims;rep(m)=1;
        kkAux=repmat(kGrid{m},rep);
        kk(:,m)=kkAux(:);
    end
    kk=reshape(kk,[NNDims NDims]);
    dimkk=ndims(kk);
    kkrad=sqrt(sum(kk.^2,dimkk));%Radial k-coordinates
    tuk=single(ones(NNDims));
    if gpu>0;tuk=gpuArray(tuk);end
    alpha=1-gibbsRing;
    if gibbsRing~=0
        fkk=0.5*(1+cos(pi*((kkrad-pi*alpha)/((1-alpha)*pi))));
        tuk(kkrad>=pi*alpha)=fkk(kkrad>=pi*alpha);
    end
    tuk(kkrad>=pi)=0;

    tuk=ifftshift(tuk);
    x=filtering(x,tuk,gpu); 
else
    NDims=min(NDims,NGR);
    for m=1:NDims
        tuk=tukeywin(N(m)+1-mod(N(m),2),gibbsRing(m))';
        tuk=tuk(1:end-(1-mod(N(m),2)));
        perm=1:NDims;perm(2)=m;perm(m)=2;
        tuk=permute(tuk,perm);
        tuk=ifftshift(tuk,m);
        x=filtering(x,tuk,gpu);
    end
end
