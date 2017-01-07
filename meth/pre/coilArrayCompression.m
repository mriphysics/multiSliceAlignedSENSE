function [SC,yC]=coilArrayCompression(S,y,perc,gpu)

%COILARRAYCOMPRESSION   Compress the measured information
%   [SC,YC] = COILARRAYCOMPRESSION(S,Y,PERC,GPU) performs coil compression
%   based on M. Buehrer, K.P. Pruessmann, P. Boesiger, and S. Kozerke.
%   ''Array Compression for MRI With Large Coil Arrays'', Magnetic
%   Resonance in Medicine, 57:1131-1139, 2007 to accelerate computations
%   S is the coil-array sensitivity map
%   Y is the measured data
%   PERC is the energy rate to be kept
%   GPU is a flag that determines whether to use gpu (1) or cpu (0)
%   computation
%   It returns:
%   SC, compressed coil-array sensitivity map
%   YC, compressed measured data
%

reg=0.001;

if perc~=1    
    %Initialize
    if ~isempty(y);NY=size(y);else NY=size(S);end
    
    %Compute P
    if gpu;S=gpuArray(S);end
    Sconj=conj(S);
    Normal=(sum(Sconj.*S,4)+reg).^(-1);
    Saux=bsxfun(@times,S,Normal);
    clear Normal
    Sconj=permute(Sconj,[1 2 3 5 4]);
    P=single(zeros(NY(4),NY(4)));
    for s=1:NY(3)
        Sbis=bsxfun(@times,Saux(:,:,s,:),Sconj(:,:,s,:,:));    
        P=P+shiftdim(sum(sum(Sbis,2),1),3);        
    end
    clear Saux Sconj Sbis
    if gpu
        P=gather(P);
    end
    
    %Compute F
    [U,F]=svd(P);
    F=diag(F); 
    
    %Thresholding
    Ftot=sum(F);
    for m=1:NY(4)
        if sum(F(1:m))/Ftot>=perc
            M=m;
            break
        end
    end
    
    %Compressing matrix
    A=ctranspose(U);
    A=permute(A(1:M,:),[2 1]);
    A=shiftdim(A,-3);
    if gpu;A=gpuArray(A);end
    
    %Compress sensitivities
    NS=size(S);    
    SC=single(zeros([NS(1:3) M]));
    if gpu
        SC=gpuArray(SC);
    end
    for m=1:M
        SC(:,:,:,m)=sum(bsxfun(@times,S,A(1,1,1,:,m)),4);
    end
    if gpu;SC=gather(SC);end
    if ~isempty(y)
        %Compress data
        yC=single(zeros([NY(1:3) M]));
        if gpu;yC=gpuArray(yC);y=gpuArray(y);end
        for m=1:M                             
            yC(:,:,:,m)=sum(bsxfun(@times,y,A(1,1,1,:,m)),4);
        end
        if gpu;yC=gather(yC);end
    else
        yC=[];
    end
else
    SC=S;yC=y;
end
