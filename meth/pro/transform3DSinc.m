function [I,IB]=transform3DSinc(I,et,di,gpu,F,FH)

%TRANSFORM3DSINC rigidly transforms volumes using sinc interpolation
%   [I,IB]=TRANSFORM3DSINC(I,ET,DI,GPU,F,FH) transforms the images with
%   sinc-based interpolation (both forwards and backwards)
%   I is the input image
%   ET are the transform factors
%   DI is a flag to indicate whether to perform direct (1) or inverse (0) 
%   transform
%   GPU is a flag that determines whether to use gpu processing
%   F is the precomputed DFT matrix (defaults to empty)
%   FH is the precomputed IDFT matrix (defaults to empty)
%   It returns:
%   I, the transformed image
%   IB, the auxiliary images after different steps


if ~exist('F','var');
    for m=1:3;F{m}=[];end
end
if ~exist('FH','var');
    for m=1:3;FH{m}=[];end
end

tr(1,:)=[1 3 2];tr(2,:)=[2 1 3];

if di==1            
    %Rotations    
    for m=1:3  
        I=fftGPU(I,tr(1,m),gpu,F{tr(1,m)});
        IB{m}=I;
        I=bsxfun(@times,et{2}{m},I);                
        I=ifftGPU(I,tr(1,m),gpu,FH{tr(1,m)});
        I=fftGPU(I,tr(2,m),gpu,F{tr(2,m)});
        I=bsxfun(@times,I,et{3}{m});
        I=ifftGPU(I,tr(2,m),gpu,FH{tr(2,m)});
        I=fftGPU(I,tr(1,m),gpu,F{tr(1,m)});
        I=bsxfun(@times,I,et{2}{m});
        if m~=3
            I=ifftGPU(I,tr(1,m),gpu,FH{tr(1,m)});
        end
    end    
    
    %Translation
    for m=1:3     
        if m~=tr(1,3)
            I=fftGPU(I,m,gpu,F{m});
        end
    end
    IB{4}=I;    
    I=bsxfun(@times,I,et{1});    
    for m=3:-1:1   
        I=ifftGPU(I,m,gpu,FH{m});       
    end
else
    %Back-translation
    for m=1:3 
        I=fftGPU(I,m,gpu,F{m});    
    end
    I=bsxfun(@times,I,et{1});%I've modified this lately    
    for m=3:-1:1  
        if m~=tr(1,3)
            I=ifftGPU(I,m,gpu,FH{m});       
        end
    end   
    %Back-rotations
    for m=3:-1:1 
        if m~=3
            I=fftGPU(I,tr(1,m),gpu,F{tr(1,m)});        
        end
        I=bsxfun(@times,I,et{2}{m});
        I=ifftGPU(I,tr(1,m),gpu,FH{tr(1,m)});
        I=fftGPU(I,tr(2,m),gpu,F{tr(2,m)});
        I=bsxfun(@times,I,et{3}{m});
        I=ifftGPU(I,tr(2,m),gpu,FH{tr(2,m)});            
        I=fftGPU(I,tr(1,m),gpu,F{tr(1,m)});
        I=bsxfun(@times,et{2}{m},I);
        if m==1
            I=sum(I,5);
        end
        I=ifftGPU(I,tr(1,m),gpu,FH{tr(1,m)});        
    end    
end
