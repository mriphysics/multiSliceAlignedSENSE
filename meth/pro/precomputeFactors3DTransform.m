function [et,etg,eth]=precomputeFactors3DTransform(xkGrid,kkGrid,kGrid,T,di,cg,gpu)

%PRECOMPUTEFACTORS3DTRANSFORM precomputes the required terms to apply a
%sinc-based rigid transform
%   [ET,ETG,ETH]=PRECOMPUTEFACTORS3DTRANSFORM(XKGRID,KKGRID,KGRID,T,DI,CG,GPU) 
%   precomputes the k-space phase multiplicative factors used for 
%   transforming the images with sinc-based interpolation.
%   XKGRID is the grid of points in the spatial-spectral domain
%   KKGRID is the grid of points in the spectral-spectral domain
%   KGRID is the grid of points in the spectral domain
%   T are the parameters of the transform
%   DI is a flag to indicate whether to perform direct (1) or inverse (0)
%   transform
%   CG is a flag to indicate the order of the derivative terms to calculate
%   GPU is a flag that determines whether to use gpu processing (defaults
%   to 0)
%   It returns:
%   ET, transform factors
%   ETG, transform gradient factors
%   ETH, transform Hessian factors
%

if ~exist('gpu','var');gpu=0;end

if gpu
    T=gpuArray(T);
    for m=1:3
        for n=1:2
            xkGrid{n}{m}=gpuArray(xkGrid{n}{m});
        end
        kGrid{m}=gpuArray(kGrid{m});        
    end
    if ~isempty(kkGrid)
        for m=1:6
            kkGrid{m}=gpuArray(kkGrid{m});
        end
    end
end

ND=ndims(T);
perm=1:ND;perm(1)=ND;perm(ND)=1;
T=permute(T,perm);
theta=T(4:6,:,:,:,:,:,:,:);
t=((-1)^di)*1i*T(1:3,:,:,:,:,:,:,:);

tantheta2=theta/2;
tantheta2=tan(tantheta2);
tantheta2j=((-1)^(di-1))*1i*tantheta2;
sintheta=sin(theta);
sintheta=((-1)^di)*1i*sintheta;    
clear T 
if cg>0
    tantheta=tan(theta);
    tanthetacuad=tantheta2.*tantheta2;
    tanthetacuad=(1+tanthetacuad)/2;
    costheta=cos(theta);    
end       
clear theta

per(1,:)=[1 3 2];per(2,:)=[2 1 3];
for m=1:3
    et{2}{m}=exp(bsxfun(@times,permute(tantheta2j(m,:,:,:,:,:,:,:),perm),xkGrid{1}{m}));%Tan exponential
    et{3}{m}=exp(bsxfun(@times,permute(sintheta(m,:,:,:,:,:,:,:),perm),xkGrid{2}{m}));%Sin exponential
end
clear tantheta2j sintheta

if cg>0  
    for m=1:3     
        etg{2}{m}=bsxfun(@times,permute(tanthetacuad(m,:,:,:,:,:,:,:),perm),1i*xkGrid{1}{m});%Tan derivative
        etg{3}{m}=bsxfun(@times,permute(costheta(m,:,:,:,:,:,:,:),perm),-1i*xkGrid{2}{m});%Sin derivative
        if cg==2
            eth{2}{m}=bsxfun(@plus,permute(tantheta2(m,:,:,:,:,:,:,:),perm),etg{2}{m});
            eth{3}{m}=bsxfun(@plus,-permute(tantheta(m,:,:,:,:,:,:,:),perm),etg{3}{m});        
        end
    end
    clear tanthetacuad costheta tantheta2 tantheta

    for m=2:3        
        for n=1:3        
            etg{m}{n}=etg{m}{n}.*et{m}{n};
            if cg==2
                eth{m}{n}=eth{m}{n}.*etg{m}{n};
            end
            etg{m}{n}=ifftshift(etg{m}{n},per(m-1,n));
            if cg==2
                eth{m}{n}=ifftshift(eth{m}{n},per(m-1,n));
            end
        end
    end
end

for m=1:3
    for n=2:3
        et{n}{m}=ifftshift(et{n}{m},per(n-1,m));
    end
end

et{1}=bsxfun(@plus,bsxfun(@plus,bsxfun(@times,permute(t(1,:,:,:,:,:,:,:),perm),kGrid{1}),bsxfun(@times,permute(t(2,:,:,:,:,:,:,:),perm),kGrid{2})),bsxfun(@times,permute(t(3,:,:,:,:,:,:,:),perm),kGrid{3}));
et{1}=exp(et{1});

clear t
if cg>0   
    for m=1:3
        etg{1}{m}=bsxfun(@times,-1i*kGrid{m},et{1});        
        for n=1:3
            etg{1}{m}=ifftshift(etg{1}{m},n);
        end
    end    
    if cg==2
        for m=1:6
            eth{1}{m}=bsxfun(@times,-kkGrid{m},et{1});
            for n=1:3
                eth{1}{m}=ifftshift(eth{1}{m},n);        
            end
        end
    end
end

for m=1:3
    et{1}=ifftshift(et{1},m);    
end
