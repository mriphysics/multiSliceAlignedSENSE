function x=resampling(x,Nres,fo,gpu,mirror)

%RESAMPLING resamples a given image using the FFT
%   X=RESAMPLING(X,NRES,FO,GPU,MIRROR) resamples an image to the number of 
%   pixels given by NRES.
%   X is the image to be resampled
%   NRES is the new resolution
%   FO determines whether the input/output is in k-space (1) or space (0,
%   default)
%   GPU determines whether to use gpu computation (1) or cpu computation
%   (0, default)
%   MIRROR -> determines whether to apply mirrored boundary conditions
%   along a given dimension to avoid ringing (defaults to no)
%   It returns:
%   X, the resampled image
%

N=size(x);
nDimsIn=length(N);
while Nres(end)==1;Nres=Nres(1:end-1);end
nDimsOu=length(Nres);

if ~exist('fo','var');fo=0;end
if ~exist('gpu','var');gpu=0;end
if ~exist('mirror','var');mirror=zeros(1,nDimsOu);end

if nDimsOu>nDimsIn
    error('Resampling dimensionality is larger than image dimensionality');
end
if nDimsIn>12
    error('Downsampling only works for arrays of up to 12 dimensions');
end
Nor=N(1:nDimsOu);
mirror=mirror(1:nDimsOu);

NorM=Nor+mirror.*Nor;
NresM=Nres+mirror.*Nres;
Nmin=min(NorM,NresM);
Nmax=max(NorM,NresM);

zeroF=ceil((Nmax+1)/2);
orig=zeroF-ceil((Nmin-1)/2);
fina=zeroF+floor((Nmin-1)/2);

for m=1:nDimsOu
    if Nor(m)~=Nres(m)
        NNres=[Nres(1:m-1) NresM(m) N(m+1:end)];        
        if m~=1
            perm=1:nDimsIn;perm(1)=m;perm(m)=1;                
            x=permute(x,perm);
            NNres=NNres(perm);
        end           
        xRes=single(zeros(NNres));
        if gpu>0
            xRes=gpuArray(xRes);
        end
        x=mirroring(x,mirror(m),1);        
        if ~fo
            x=fftGPU(x,1,gpu)/NorM(m);
        end
        x=fftshift(x,1);
        if Nor(m)>Nres(m)
            xRes=x(orig(m):fina(m),:,:,:,:,:,:,:,:,:,:,:);
        else
            xRes(orig(m):fina(m),:,:,:,:,:,:,:,:,:,:,:)=x;
        end       
        x=xRes;
        x=ifftshift(x,1);
        if ~fo
            x=ifftGPU(x,1,gpu)*NresM(m);
        end
        x=mirroring(x,mirror(m),0);
        if m~=1
            x=permute(x,perm);
        end 
    end
end
