function x=ifftGPU(x,m,gpu,FH)

%IFFTGPU   Optimal arrangement of GPU-based IFFT computation
%   X=IFFTGPU(X,M,GPU,FH) applies IFFT's over rearranged arrays so that the 
%   IDFT is applied across the first dimension. This has accelerated GPU 
%   computations in our setting
%   * X is the array on which to apply the IDFT
%   * M is the direction across which to apply the IDFT
%   * GPU is a flag that determines whether to use matrix-based gpu (2), 
%   gpu (1) or cpu (0) computation
%   * F is the IDFT matrix
%   It returns:
%   * X, the IDFT-transformed array
%

if ~exist('FH','var');FH=[];end

if gpu==0
    x=ifft(x,[],m);
elseif gpu==1
    if m==1
        if size(x,1)~=1
            x=ifft(x);
        end
    else
        perm=1:ndims(x);perm(1)=m;perm(m)=1;
        x=permute(x,perm);
        if size(x,1)~=1
            x=ifft(x);
        end
        x=permute(x,perm);
    end
else
    N=size(x,m);
    if N~=1
        if isempty(FH)
            FH=conj(gpuArray(single(dftmtx(N))))/N;  
        end     
        if m~=1
            perm=1:ndims(x);perm(1)=m;perm(m)=1;        
            x=permute(x,perm);
        end
        if ~ismatrix(x)
            S=size(x);S(end+1:2)=1;        
            x=reshape(x,[S(1) prod(S(2:end))]);
        end
        x=FH*x;
        if exist('S','var')
            x=reshape(x,S);
        end
        if m~=1
            x=permute(x,perm);     
        end
    end
end