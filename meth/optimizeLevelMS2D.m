function [x,T,outlD]=optimizeLevelMS2D(x,y,T,S,W,Ak,xGrid,mNorm,parX,debug,res,estT,gpu,SlTh,SlOv,outlD,alpha)

%OPTIMIZELEVELMS2D gets the aligned reconstruction for a given resolution
%level
%   [X,T,OUTLD]=OPTIMIZELEVELMS2D(X,Y,T,S,W,AK,XGRID,MNORM,PARX,DEBUG,RES,ESTT,GPU,SLTH,SLOV,OUTLD,ALPHA)
%   applies the motion corrected reconstruction algorithm for MS sequences
%   at a given resolution level.
%   X is the input image
%   Y is the measured data
%   T are the transform parameters
%   S are the sensitivities
%   W is a spatial mask
%   AK is a sampling mask in the phase encoding direction
%   XGRID is a grid of points in the spatial domain
%   MNORM is a normalization factor for measurements
%   PARX is a structure with the reconstruction parameters
%   DEBUG is a flag to print debug information: 0-> no; 1-> basic; 2-> 
%   verbose
%   RES is the subsampling factor for this resolution level
%   ESTT is a flag that determines whether to apply motion correction
%   GPU is a flag that determines whether to use gpu processing
%   SLTH is the slice thickness
%   SLOV is the slice overlap
%   OUTLD is a mask for shot rejection
%   ALPHA is the regularization factor to treat the slice profiles
%   It returns:
%   X, the reconstructed image
%   T, the estimated tranform paramaters
%   OUTLD, the estimated shot rejection mask

%Parameters fixed in this high level implementation
reg=0.001;

%ROI computation
ROIsec=0;
ROI=computeROI(W,ROIsec);
if parX.threeD
    if parX.correct==0;NzSlab=3;else NzSlab=7;end
else
    NzSlab=1;
end

%z-Overlap
H=sliceProfile(NzSlab,SlTh,SlOv,parX);

Ha=zeros(1,size(W,3));Hb=zeros(1,size(W,3));
if ~parX.threeD;alpha=0;end
H(1)=1;
Ha(1)=1;Ha(2)=-1;
Hb(1)=1;Hb(end)=-1;
Ha=fft(Ha);
Hb=fft(Hb);
Hsmooth=alpha*(conj(Ha).*Ha+conj(Hb).*Hb);%This is what is denoted W in the paper
Hsmooth=permute(Hsmooth,[1 3 2]);

%ROI extraction in the readout direction
y=extractROI(y,ROI,2,1);S=extractROI(S,ROI,2,1);W=extractROI(W,ROI,2,1);x=extractROI(x,ROI,2,1);xGrid{2}=extractROI(xGrid{2},ROI,2,1);

N=size(W);
N(3)=NzSlab;    

kGrid{1}=permute(single(-floor(N(1)/2):ceil(N(1)/2)-1),[2 1])/res;
kGrid{2}=single(-floor(N(2)/2):ceil(N(2)/2)-1)/res;
kGrid{3}=permute(single(-floor(N(3)/2):ceil(N(3)/2)-1),[1 3 2]);

%%Update transform terms
for m=1:3
    kGrid{m}=2*pi*kGrid{m}/N(m);
end    
yZ{2}{2}=single(zeros([1 1 N(3) 1 1 size(xGrid{3},3)]));
xGrid{3}=extractSlabs(xGrid{3},N(3),1,1,yZ);

per(1,:)=[1 3 2];per(2,:)=[2 1 3];
xkGrid=cell(2,3);
for n=1:2
    for m=1:3
        xkGrid{n}{m}=bsxfun(@times,xGrid{per(3-n,m)},kGrid{per(n,m)});%First index denotes that the k takes the first dimension of xkgrid; second index denotes the rotation        
    end
end
fact(1,:)=[1 2 3 1 1 2];fact(2,:)=[1 2 3 2 3 3];
kkGrid=cell(1,6);
for m=1:6
    kkGrid{m}=bsxfun(@times,kGrid{fact(1,m)},kGrid{fact(2,m)});
end

%Initialize solver
y=y/mNorm;x=x/mNorm;
winic=1;
NT=size(T);
w=winic*ones(NT(1:6));
flagw=zeros(NT(1:6));
nExtern=1000;
nX=5;
nT=1;

%Preconditioner computation
SH=conj(S);
Precond=sum(real(SH.*S),4);
Precond=(Precond+reg).^(-1);

y=fft(y);
for n=1:nExtern       
    if debug==2;tstart=tic;end
    xant=x;
    %Solve for x
    %Block sizes for gpu computation
    if res==1 && sum(sum(sum(T)))~=0    
        BlSz(1)=1;BlSz(2)=ceil(size(T,6)/2);
    else
        BlSz(1)=floor(size(T,5)/2);BlSz(2)=size(T,6);
        if BlSz(1)==0;BlSz(1)=1;end
    end    
    if parX.correct>0 && estT
        x=solveXMS2D(x,y,W,T,H,Hsmooth,S,SH,Precond,Ak,xkGrid,kGrid,nX,0,NzSlab,gpu,outlD,BlSz);
    else
        x=solveXMS2D(x,y,W,T,H,Hsmooth,S,SH,Precond,Ak,xkGrid,kGrid,inf,parX.toler,NzSlab,gpu,outlD,BlSz);
    end
    if debug==2;telapsed=toc(tstart);fprintf('Time solving x: %.4fs\n',telapsed);end

    %Solve for T    
    if parX.correct>0 && estT
        if debug==2;tstart=tic;end
        BlSz=1;
        [T,w,flagw,outlD]=solveTMS2D(x,y,T,H,S,Ak,xkGrid,kkGrid,kGrid,nT,w,flagw,NzSlab,gpu,parX.outlP,parX.thplc,BlSz);
        if debug==2;telapsed=toc(tstart);fprintf('Time solving T: %.4fs\n',telapsed);end
        
        xant=x-xant;
        xant=real(xant.*conj(xant));
        xant=max(xant(:));       
        if xant<parX.toler
            if debug>0;fprintf('Iteration XT %04d - Error %0.2g \n',n,xant);end 
            break
        end                     
    else
        break
    end
end

x=x*mNorm;
x=extractROI(x,ROI,2,0);
