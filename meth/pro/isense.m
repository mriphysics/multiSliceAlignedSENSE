function x=isense(x,m,NS,NY,FOV)

%ISENSE   Inverse SENSE folding operation
%   X=ISENSE(X,M,NS,NY,FOV) applies the inverse SENSE folding operation
%   along a given dimension of the input array
%   * X is the array on which to operate
%   * M is the direction across which to operate
%   * NS is the spatial size of the array
%   * NY is the spectral size of the array
%   * FOV indicates the size of the slab to be folded
%   It returns:
%   * X, the inverse SENSE folded array
%

oddFactSENSE=2*ceil((ceil(NS/NY)-1)/2)+1;
oFRed=oddFactSENSE-2;

nDimsIn=ndims(x);

if NS~=NY
    %if mod(NS,NY)~=0
        if m~=1
            perm=1:nDimsIn;perm(1)=m;perm(m)=1;
            x=permute(x,perm);
        end          
        NX=size(x);NX(end+1:2);
        x=reshape(x,[NX(1) prod(NX(2:end))]);
        indUnf=[1+FOV(1):NX(1) repmat(1:NX(1),[1 oFRed]) 1:NX(1)-FOV(2)];
        x=x(indUnf,:);
        x=reshape(x,[NS NX(2:end)]);   
        if m~=1
            x=permute(x,perm);
        end
    %else
    %    perm=zeros(1,nDimsIn);perm(m)=-FOV(1);        
    %    x=circshift(x,perm);
    %    perm=ones(1,nDimsIn);perm(m)=NS/NY;
    %    x=repmat(x,perm);                
    %end
end
