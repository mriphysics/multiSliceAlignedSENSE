function [G,GB,GC]=transform3DSincGradient(xB,et,etg,full,gpu,F,FH)

%TRANSFORM3DSINCGRADIENT computes the gradient of the rigid transform
%   [G,GB,GC]=TRANSFORM3DSINCGRADIENT(xB,ET,ETG,FULL,GPU,F,FH) obtains the
%   gradient of the transform of the images
%   xB are the images before the first, second and third rotations and
%   before the translation
%   ET are the transform factors
%   ETG are the transform gradient factors
%   FULL indicates whether the gradient along different transform
%   parameters is stored as an array (1) or as a cell (0)
%   GPU is a flag that determines whether to use gpu processing
%   F is the precomputed DFT matrix (defaults to empty)
%   FH is the precomputed IDFT matrix (defaults to empty)
%   It returns:
%   G, the gradient of the transformed image
%   GB, the gradient of the first, second and third rotations before 
%   applying the translation
%   GC, the gradient of the first, first and second rotations before
%   applying the second, third and third rotations respectively
%

if ~exist('F','var');
    for m=1:3;F{m}=[];end
end
if ~exist('FH','var');
    for m=1:3;FH{m}=[];end
end

GB=[];GC=[];

%Translation parameters
for m=1:3    
    x{1}=bsxfun(@times,xB{4},etg{1}{m});
    for n=1:3       
        x{1}=ifftGPU(x{1},n,gpu,FH{n});
    end    
    G{m}=x{1};
end

%First rotation
x{1}=bsxfun(@times,xB{1},et{2}{1});
x{2}=bsxfun(@times,xB{1},etg{2}{1});
for m=1:2
    x{m}=ifftGPU(x{m},1,gpu,FH{1});
    x{m}=fftGPU(x{m},2,gpu,F{2});
end
x{2}=bsxfun(@times,etg{3}{1},x{1})+bsxfun(@times,et{3}{1},x{2});
x{1}=bsxfun(@times,et{3}{1},x{1});    
for m=1:2
    x{m}=ifftGPU(x{m},2,gpu,FH{2});       
    x{m}=fftGPU(x{m},1,gpu,F{1});
end
x{1}=bsxfun(@times,etg{2}{1},x{1})+bsxfun(@times,et{2}{1},x{2});
x{1}=ifftGPU(x{1},1,gpu,FH{1});

x{1}=fftGPU(x{1},3,gpu,F{3});
GC{1}=x{1};
x{1}=bsxfun(@times,x{1},et{2}{2});
x{1}=ifftGPU(x{1},3,gpu,FH{3});
x{1}=fftGPU(x{1},1,gpu,F{1});
x{1}=bsxfun(@times,x{1},et{3}{2});
x{1}=ifftGPU(x{1},1,gpu,FH{1});
x{1}=fftGPU(x{1},3,gpu,F{3});
x{1}=bsxfun(@times,x{1},et{2}{2});
x{1}=ifftGPU(x{1},3,gpu,FH{3});

x{1}=fftGPU(x{1},2,gpu,F{2});
GC{2}=x{1};
x{1}=bsxfun(@times,x{1},et{2}{3});
x{1}=ifftGPU(x{1},2,gpu,FH{2});
x{1}=fftGPU(x{1},3,gpu,F{3});
x{1}=bsxfun(@times,x{1},et{3}{3});
x{1}=ifftGPU(x{1},3,gpu,FH{3});
x{1}=fftGPU(x{1},2,gpu,F{2});
x{1}=bsxfun(@times,x{1},et{2}{3});

for m=1:2:3
    x{1}=fftGPU(x{1},m,gpu,F{m});
end
GB{1}=x{1};
x{1}=bsxfun(@times,x{1},et{1});
for m=3:-1:1    
    x{1}=ifftGPU(x{1},m,gpu,FH{m});   
end
G{4}=x{1};


%Second rotation
x{1}=bsxfun(@times,xB{2},et{2}{2});
x{2}=bsxfun(@times,xB{2},etg{2}{2});
for m=1:2
    x{m}=ifftGPU(x{m},3,gpu,FH{3});
    x{m}=fftGPU(x{m},1,gpu,F{1});
end
x{2}=bsxfun(@times,etg{3}{2},x{1})+bsxfun(@times,et{3}{2},x{2});
x{1}=bsxfun(@times,et{3}{2},x{1});
for m=1:2
    x{m}=ifftGPU(x{m},1,gpu,FH{1});
    x{m}=fftGPU(x{m},3,gpu,F{3});
end
x{1}=bsxfun(@times,etg{2}{2},x{1})+bsxfun(@times,et{2}{2},x{2});
x{1}=ifftGPU(x{1},3,gpu,FH{3});

x{1}=fftGPU(x{1},2,gpu,F{2});
GC{3}=x{1};
x{1}=bsxfun(@times,x{1},et{2}{3});
x{1}=ifftGPU(x{1},2,gpu,FH{2});
x{1}=fftGPU(x{1},3,gpu,F{3});
x{1}=bsxfun(@times,x{1},et{3}{3});
x{1}=ifftGPU(x{1},3,gpu,FH{3});
x{1}=fftGPU(x{1},2,gpu,F{2});
x{1}=bsxfun(@times,x{1},et{2}{3});

for m=1:2:3
    x{1}=fftGPU(x{1},m,gpu,F{m});
end
GB{2}=x{1};
x{1}=bsxfun(@times,x{1},et{1});
for m=3:-1:1    
    x{1}=ifftGPU(x{1},m,gpu,FH{m});    
end
G{5}=x{1};


%Third rotation
x{1}=bsxfun(@times,xB{3},et{2}{3});
x{2}=bsxfun(@times,xB{3},etg{2}{3});
for m=1:2
    x{m}=ifftGPU(x{m},2,gpu,FH{2});
    x{m}=fftGPU(x{m},3,gpu,F{3});
end
x{2}=bsxfun(@times,etg{3}{3},x{1})+bsxfun(@times,et{3}{3},x{2});
x{1}=bsxfun(@times,et{3}{3},x{1});
for m=1:2
    x{m}=ifftGPU(x{m},3,gpu,FH{3});
    x{m}=fftGPU(x{m},2,gpu,F{2});
end
x{1}=bsxfun(@times,etg{2}{3},x{1})+bsxfun(@times,et{2}{3},x{2});

for m=1:2:3
    x{1}=fftGPU(x{1},m,gpu,F{m});        
end
GB{3}=x{1};
x{1}=bsxfun(@times,x{1},et{1});
for m=3:-1:1       
    x{1}=ifftGPU(x{1},m,gpu,FH{m});        
end
G{6}=x{1};

if full
    G=cat(6,G{1},G{2},G{3},G{4},G{5},G{6});
end
