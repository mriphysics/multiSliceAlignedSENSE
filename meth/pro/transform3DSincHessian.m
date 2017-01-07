function G=transform3DSincHessian(xB,GB,GC,et,etg,eth,mH,gpu,F,FH)

%TRANSFORM3DSINCHESSIAN computes the Hessian of the rigid transform
%   G=TRANSFORM3DSINCHESSIAN(XB,GB,GC,ET,ETG,ETH,MH,GPU,F,FH) obtains the
%   Hessian of the transform of the images
%   xB are the images before the first, second and third rotations and
%   before the translation
%   GB is the gradient of the first, second and third rotations before 
%   applying the translation
%   GC is the gradient of the first, first and second rotations before
%   applying the second, third and third rotations respectively
%   ET are the transform factors
%   ETG are the transform gradient factors
%   ETH are the transform Hessian factors
%   MH is the component of the Hessian
%   GPU is a flag that determines whether to use gpu processing
%   F is the precomputed DFT matrix (defaults to empty)
%   FH is the precomputed IDFT matrix (defaults to empty)
%   It returns:
%   G, the Hessian of the transformed image
%

if ~exist('F','var');
    for m=1:3;F{m}=[];end
end
if ~exist('FH','var');
    for m=1:3;FH{m}=[];end
end

%Translation parameters
for m=1:6
    if m==mH    
        x{1}=bsxfun(@times,xB{4},eth{1}{m});%I've modified this lately
        for n=1:3            
            x{1}=ifftGPU(x{1},n,gpu,FH{n});            
        end  
        G=x{1};
    end
end


%Translation-rotation parameters
%First rotation
for n=1:3
    for m=1:3
        if (6+3*(n-1)+m)==mH
            x{1}=bsxfun(@times,GB{n},etg{1}{m});
            for o=1:3              
                x{1}=ifftGPU(x{1},o,gpu,FH{o});                
            end   
            G=x{1};
        end
    end 
end

%Rotation cross terms
%First-second
if mH==16
    x{1}=bsxfun(@times,GC{1},et{2}{2});
    x{2}=bsxfun(@times,GC{1},etg{2}{2});
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
    x{1}=bsxfun(@times,x{1},et{1});
    for m=3:-1:1     
        x{1}=ifftGPU(x{1},m,gpu,FH{m});        
    end
    G=x{1};
end


%First-third
if mH==17
    x{2}=bsxfun(@times,GC{2},etg{2}{3});
    x{1}=bsxfun(@times,GC{2},et{2}{3});
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
    x{1}=bsxfun(@times,x{1},et{1});
    for m=3:-1:1        
        x{1}=ifftGPU(x{1},m,gpu,FH{m});        
    end
    G=x{1};
end

%Second-third
if mH==18
    x{2}=bsxfun(@times,GC{3},etg{2}{3});
    x{1}=bsxfun(@times,GC{3},et{2}{3});
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
    x{1}=bsxfun(@times,x{1},et{1});
    for m=3:-1:1      
        x{1}=ifftGPU(x{1},m,gpu,FH{m});       
    end
    G=x{1};
end


%Rotation second order
%First rotation
if mH==19
    x{1}=bsxfun(@times,xB{1},et{2}{1});
    x{2}=bsxfun(@times,xB{1},eth{2}{1});
    x{3}=bsxfun(@times,xB{1},etg{2}{1});
    x{4}=bsxfun(@times,xB{1},etg{2}{1});
    for m=1:4
        x{m}=ifftGPU(x{m},1,gpu,FH{1});
        x{m}=fftGPU(x{m},2,gpu,F{2});
    end
    x{2}=bsxfun(@times,eth{3}{1},x{1})+bsxfun(@times,et{3}{1},x{2});
    x{5}=bsxfun(@times,etg{3}{1},x{1})+bsxfun(@times,et{3}{1},x{3});
    x{6}=bsxfun(@times,etg{3}{1},x{1})+bsxfun(@times,et{3}{1},x{4});    
    x{1}=bsxfun(@times,et{3}{1},x{1});
    x{3}=bsxfun(@times,etg{3}{1},x{3});
    x{4}=bsxfun(@times,etg{3}{1},x{4});
    for m=1:6
        x{m}=ifftGPU(x{m},2,gpu,FH{2});
        x{m}=fftGPU(x{m},1,gpu,F{1});
    end
    x{1}=bsxfun(@times,eth{2}{1},x{1})+bsxfun(@times,et{2}{1},x{2});
    x{2}=bsxfun(@times,etg{2}{1},x{5})+bsxfun(@times,et{2}{1},x{3});
    x{3}=bsxfun(@times,etg{2}{1},x{6})+bsxfun(@times,et{2}{1},x{4});
    x{1}=x{1}+x{2}+x{3};
    x{1}=ifftGPU(x{1},1,gpu,FH{1});
        
    x{1}=fftGPU(x{1},3,gpu,F{3});
    x{1}=bsxfun(@times,x{1},et{2}{2});
    x{1}=ifftGPU(x{1},3,gpu,FH{3});
    x{1}=fftGPU(x{1},1,gpu,F{1});
    x{1}=bsxfun(@times,x{1},et{3}{2});
    x{1}=ifftGPU(x{1},1,gpu,FH{1});
    x{1}=fftGPU(x{1},3,gpu,F{3});
    x{1}=bsxfun(@times,x{1},et{2}{2});
    x{1}=ifftGPU(x{1},3,gpu,FH{3});
    
    x{1}=fftGPU(x{1},2,gpu,F{2});
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
    x{1}=bsxfun(@times,x{1},et{1});
    for m=3:-1:1      
        x{1}=ifftGPU(x{1},m,gpu,FH{m});        
    end
    G=x{1};
end


%Second rotation
if mH==20
    x{1}=bsxfun(@times,xB{2},et{2}{2});
    x{2}=bsxfun(@times,xB{2},eth{2}{2});
    x{3}=bsxfun(@times,xB{2},etg{2}{2});
    x{4}=bsxfun(@times,xB{2},etg{2}{2});
    for m=1:4
        x{m}=ifftGPU(x{m},3,gpu,FH{3});
        x{m}=fftGPU(x{m},1,gpu,F{1});
    end
    x{2}=bsxfun(@times,eth{3}{2},x{1})+bsxfun(@times,et{3}{2},x{2});
    x{5}=bsxfun(@times,etg{3}{2},x{1})+bsxfun(@times,et{3}{2},x{3});
    x{6}=bsxfun(@times,etg{3}{2},x{1})+bsxfun(@times,et{3}{2},x{4});    
    x{1}=bsxfun(@times,et{3}{2},x{1});
    x{3}=bsxfun(@times,etg{3}{2},x{3});
    x{4}=bsxfun(@times,etg{3}{2},x{4});
    for m=1:6
        x{m}=ifftGPU(x{m},1,gpu,FH{1});
        x{m}=fftGPU(x{m},3,gpu,F{3});
    end
    x{1}=bsxfun(@times,eth{2}{2},x{1})+bsxfun(@times,et{2}{2},x{2});
    x{2}=bsxfun(@times,etg{2}{2},x{5})+bsxfun(@times,et{2}{2},x{3});
    x{3}=bsxfun(@times,etg{2}{2},x{6})+bsxfun(@times,et{2}{2},x{4});
    x{1}=x{1}+x{2}+x{3};
    x{1}=ifftGPU(x{1},3,gpu,FH{3});
        
    x{1}=fftGPU(x{1},2,gpu,F{2});
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
    x{1}=bsxfun(@times,x{1},et{1});
    for m=3:-1:1       
        x{1}=ifftGPU(x{1},m,gpu,FH{m});        
    end
    G=x{1};
end


%Third rotation
if mH==21
    x{1}=bsxfun(@times,xB{3},et{2}{3});
    x{2}=bsxfun(@times,xB{3},eth{2}{3});
    x{3}=bsxfun(@times,xB{3},etg{2}{3});
    x{4}=bsxfun(@times,xB{3},etg{2}{3});
    for m=1:4
        x{m}=ifftGPU(x{m},2,gpu,FH{2});
        x{m}=fftGPU(x{m},3,gpu,F{3});
    end
    x{2}=bsxfun(@times,eth{3}{3},x{1})+bsxfun(@times,et{3}{3},x{2});
    x{5}=bsxfun(@times,etg{3}{3},x{1})+bsxfun(@times,et{3}{3},x{3});
    x{6}=bsxfun(@times,etg{3}{3},x{1})+bsxfun(@times,et{3}{3},x{4});    
    x{1}=bsxfun(@times,et{3}{3},x{1});
    x{3}=bsxfun(@times,etg{3}{3},x{3});
    x{4}=bsxfun(@times,etg{3}{3},x{4});
    for m=1:6
        x{m}=ifftGPU(x{m},3,gpu,FH{3});
        x{m}=fftGPU(x{m},2,gpu,F{2});
    end
    x{1}=bsxfun(@times,eth{2}{3},x{1})+bsxfun(@times,et{2}{3},x{2});
    x{2}=bsxfun(@times,etg{2}{3},x{5})+bsxfun(@times,et{2}{3},x{3});
    x{3}=bsxfun(@times,etg{2}{3},x{6})+bsxfun(@times,et{2}{3},x{4});
    x{1}=x{1}+x{2}+x{3};

    for m=1:2:3
        x{1}=fftGPU(x{1},m,gpu,F{m});        
    end
    x{1}=bsxfun(@times,x{1},et{1});
    for m=3:-1:1        
        x{1}=ifftGPU(x{1},m,gpu,FH{m});        
    end
    G=x{1};
end
