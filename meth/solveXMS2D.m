function x=solveXMS2D(x,y,W,T,H,Hsmooth,S,SH,Precond,Ak,xkGrid,kGrid,nX,toler,NzSlab,gpu,outlD,BlSz)

%SOLVEXMS2D reconstructs an image under motion
%   X=SOLVEXMS2D(X,Y,W,T,H,HSMOOTH,S,SH,PRECOND,AK,XKGRID,KGRID,NX,TOLER,NZSLAB,GPU,OUTLD,BLSZ)
%   computes the best x for a given T.
%   Inputs:
%   X is the input image
%   Y is the measured data
%   W is a spatial mask
%   T are the transform parameters
%   H is the slice profile filter
%   HSMOOTH is the regularization kernel
%   S are the sensitivities
%   SH are the complex conjugate of the sensitivities
%   PRECOND is the preconditioner
%   AK is a sampling mask in the phase encoding direction
%   XKGRID is a grid of points in the spatio-spectral domain
%   KGRID is a grid of points in the spectral domain
%   NX is the maximum number of iterations of CG method
%   TOLER is the relative error for convergence of CG method
%   NZSLAB is the number of slices of the slab
%   GPU is a flag that determines whether to use gpu processing
%   OUTLD is a mask for shot rejection
%   BLSZ is the size of blocks for gpu processing
%   It returns,
%   X, the reconstructed image
%

if gpu==1;gpuF=2;else gpuF=0;end

N=size(x);NT=size(T);
for p=1:2
    NRun(p)=ceil(NT(4+p)/BlSz(p));
    NRunRem(p)=mod(NT(4+p),BlSz(p));
    for s=1:NRun(p)
        if s~=NRun(p) || NRunRem(p)==0
            vS{p}{s}=(s-1)*BlSz(p)+1:s*BlSz(p);
        else
            vS{p}{s}=(s-1)*BlSz(p)+1:(s-1)*BlSz(p)+NRunRem(p);
        end
    end
end

if sum(T(:))~=0
    etDir=cell(NRun(1),NRun(2));etInv=cell(NRun(1),NRun(2));
    for s=1:NRun(1)
        for t=1:NRun(2)
            xkGridAux=xkGrid;
            xkGridAux{1}{3}=xkGridAux{1}{3}(:,:,:,:,:,vS{2}{t});
            xkGridAux{2}{2}=xkGridAux{2}{2}(:,:,:,:,:,vS{2}{t});
            etDir{s}{t}=precomputeFactors3DTransform(xkGridAux,[],kGrid,T(:,:,:,:,vS{1}{s},vS{2}{t},:),1,0);
            etInv{s}{t}=precomputeFactors3DTransform(xkGridAux,[],kGrid,T(:,:,:,:,vS{1}{s},vS{2}{t},:),0,0);             
        end
    end
end

%CG method
%Initialization
NS=size(S);NY=size(y);
over=(NS(1)-NY(1))/2;
disc=(3*NY(1)-NS(1))/2;
FOV=[floor(over) ceil(over)];iFOV=[floor(disc) ceil(disc)];

if gpu;Ak=gpuArray(Ak);y=gpuArray(y);S=gpuArray(S);SH=gpuArray(SH);W=gpuArray(W);outlD=gpuArray(outlD);x=gpuArray(x);end
AkOutlDisc=bsxfun(@times,Ak,outlD);

if gpu
    yEnd=gpuArray(single(zeros([N(1:2) NzSlab 1 1 N(3)]))); 
else
    yEnd=single(zeros([N(1:2) NzSlab 1 1 N(3)])); 
end

for s=1:NRun(1)
    if sum(T(:))~=0
        yS=bsxfun(@times,y,outlD(:,:,:,:,vS{1}{s}));
        yS=bsxfun(@times,yS,Ak(:,:,:,:,vS{1}{s}));
    else
        yS=y;
    end    
    yS=ifftGPU(yS,1,gpuF);
    yS=isense(yS,1,NS(1),NY(1),iFOV);
    yS=sum(bsxfun(@times,yS,SH),4); 
    NZR=size(yS);NZR(end+1:5)=1;
    if s==1
        yZR{1}{2}=single(zeros(NZR));
        yZR{2}{1}=single(zeros([NZR(1:2) NzSlab NZR(4:5) NZR(3)]));
        if gpu;yZR{2}{1}=gpuArray(yZR{2}{1});yZR{1}{2}=gpuArray(yZR{1}{2});end
    end
    yS=extractSlabs(yS,NzSlab,1,0,yZR);
    yS=filtering(yS,H,gpuF);
    if s==1
        NX=size(yS);    
        for m=1:3
            F{m}=gpuArray(single(dftmtx(NX(m))));
            FH{m}=conj(F{m})/NX(m);  
        end
    end
    
    for t=1:NRun(2)
        if sum(T(:))~=0
            gpuTSt(etInv);              
            yEnd(:,:,:,:,:,vS{2}{t})=yEnd(:,:,:,:,:,vS{2}{t})+transform3DSinc(yS(:,:,:,:,:,vS{2}{t}),etS,0,gpuF,F,FH);
        else
            yEnd(:,:,:,:,:,vS{2}{t})=yEnd(:,:,:,:,:,vS{2}{t})+sum(yS(:,:,:,:,:,vS{2}{t}),5);
        end
    end    
end
y=yEnd;
clear yEnd yS
NZR=size(y);
NZR(end+1:6)=1;
yZR{1}{1}=single(zeros([NZR(1:2) NZR(6)]));
yZR{2}{2}=single(zeros(NZR));
if gpu;yZR{1}{1}=gpuArray(yZR{1}{1});yZR{2}{2}=gpuArray(yZR{2}{2});end
y=extractSlabs(y,NzSlab,0,0,yZR);
y=W.*y;

Ap=applyCG(x);
r=y-Ap; 
z=Precond.*r;
p=z; 
rsold=sum(sum(sum(conj(z).*r)));

%Iterations
n=1;
while 1
    Ap=applyCG(p);
    al=conj(rsold)/sum(sum(sum(conj(p).*Ap)));
    xup=al*p;
    x=x+xup;
    xup=real(xup.*conj(xup));
    xup=max(max(max(xup)));
    if xup<toler || n>=nX
        if toler~=0;fprintf('Iteration CG %04d - Error %0.2g \n',n,xup);end
        break
    end
    r=r-al*Ap;
    z=Precond.*r;
    rsnew=sum(sum(sum(conj(z).*r)));
    be=rsnew/rsold;
    p=z+be*p;
    rsold=rsnew;  
    if sqrt(abs(rsnew))<1e-10
        break;
    end
    n=n+1;
end

if gpu;x=gather(x);end

function x=applyCG(x) 
    xB=filtering(x,Hsmooth,gpuF);
    x=extractSlabs(x,NzSlab,1,1,yZR);
    NX=size(x);
    xEnd=single(zeros(NX));
    if gpu;xEnd=gpuArray(xEnd);end
    for s=1:NRun(1)
        if sum(T(:))~=0
            xS=single(zeros([NX(1:4) length(vS{1}{s}) NX(6)]));
        else
            xS=single(zeros(NX));
        end
        if gpu;xS=gpuArray(xS);end
        for t=1:NRun(2)
            if sum(T(:))~=0    
                gpuTSt(etDir);
                xS(:,:,:,:,:,vS{2}{t})=transform3DSinc(x(:,:,:,:,:,vS{2}{t}),etS,1,gpuF,F,FH);
            else
                xS(:,:,:,:,:,vS{2}{t})=x(:,:,:,:,:,vS{2}{t});
            end                
        end
        xS=filtering(xS,H,gpuF);        
        xS=extractSlabs(xS,NzSlab,0,1,yZR);
        xS=bsxfun(@times,xS,S);
        xS=sense(xS,1,NS(1),NY(1),FOV);
        if sum(T(:))~=0
            xS=fftGPU(xS,1,gpuF);    
            xS=bsxfun(@times,xS,AkOutlDisc(:,:,:,:,vS{1}{s}));   
            xS=ifftGPU(xS,1,gpuF);
        end
        xS=isense(xS,1,NS(1),NY(1),iFOV);                     
        xS=sum(bsxfun(@times,xS,SH),4);
        xS=extractSlabs(xS,NzSlab,1,0,yZR);
        xS=filtering(xS,H,gpuF);
        for t=1:NRun(2)
            if sum(T(:))~=0 
                gpuTSt(etInv);               
                xEnd(:,:,:,:,:,vS{2}{t})=xEnd(:,:,:,:,:,vS{2}{t})+transform3DSinc(xS(:,:,:,:,:,vS{2}{t}),etS,0,gpuF,F,FH); 
            else
                xEnd(:,:,:,:,:,vS{2}{t})=xEnd(:,:,:,:,:,vS{2}{t})+sum(xS(:,:,:,:,:,vS{2}{t}),5);
            end
        end        
    end
    x=xEnd;
    x=extractSlabs(x,NzSlab,0,0,yZR);
    x=x+xB;
    x=x.*W;
end

function gpuTSt(etIn)
    if gpu
        etS{1}=gpuArray(etIn{s}{t}{1});
        for m=1:3
            for l=2:3
                etS{l}{m}=gpuArray(etIn{s}{t}{l}{m});            
            end
        end
    else
        etS{1}=etIn{s}{t}{1};
        for m=1:3
            for l=2:3
                etS{l}{m}=etIn{s}{t}{l}{m};
            end
        end
    end
end

end