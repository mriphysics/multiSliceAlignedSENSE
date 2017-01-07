function x=fftGPU(x,m,gpu,F)

%FFTGPU   Optimal arrangement of GPU-based FFT computation
%   X=FFTGPU(X,M,GPU,F) applies FFT's over rearranged arrays so that the 
%   DFT is applied across the first dimension. This has accelerated GPU 
%   computations in our setting
%   * X is the array on which to apply the DFT
%   * M is the direction across which to apply the DFT
%   * GPU is a flag that determines whether to use matrix-based gpu (2), 
%   gpu (1) or cpu (0) computation
%   * F is the DFT matrix
%   It returns:
%   * X, the DFT-transformed array
%

if ~exist('F','var');F=[];end

if gpu==0
    x=fft(x,[],m);
elseif gpu==1
    if m==1
        if size(x,1)~=1
            x=fft(x);
        end
    else
        perm=1:ndims(x);perm(1)=m;perm(m)=1;
        x=permute(x,perm);
        if size(x,1)~=1
            x=fft(x);
        end
        x=permute(x,perm);
    end
else
    N=size(x,m);
    if N~=1
        if isempty(F)
            F=gpuArray(single(dftmtx(N)));            
        end
        if m~=1
            perm=1:ndims(x);perm(1)=m;perm(m)=1;
            x=permute(x,perm);
        end
        if ~ismatrix(x)
            S=size(x);S(end+1:2)=1;
            x=reshape(x,[S(1) prod(S(2:end))]);
        end
        x=F*x;
        if exist('S','var')
            x=reshape(x,S);
        end
        if m~=1
            x=permute(x,perm);            
        end
    end
end