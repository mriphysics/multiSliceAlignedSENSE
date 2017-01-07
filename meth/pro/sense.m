function xb=sense(x,m,NS,NY,FOV)

%ISENSE   SENSE folding operation
%   X=SENSE(X,M,NS,NY,FOV) applies the SENSE folding operation along a 
%   given dimension of the input array
%   * X is the array on which to operate
%   * M is the direction across which to operate
%   * NS is the spatial size of the array
%   * NY is the spectral size of the array
%   * FOV indicates the size of the slab to be folded
%   It returns:
%   * X, the inverse SENSE folded array
%

oddFactSENSE=2*ceil((ceil(NS/NY)-1)/2)+1;
oFRed=(oddFactSENSE-3)/2;
oFRedT=oFRed*NY;

nDimsIn=ndims(x);

if NS~=NY   
    %if mod(NS,NY)~=0%Slower                         
         if m~=1
             perm=1:nDimsIn;perm(1)=m;perm(m)=1;
             x=permute(x,perm);
         end
         N=size(x);
         x=reshape(x,[N(1) prod(N(2:end))]);        
         xb=x(FOV(2)+1:NS-FOV(1),:);
         for s=1:oFRed
             xb=xb+x(FOV(2)+1-s*NY:NS-FOV(1)-s*NY,:)+x(FOV(2)+1+s*NY:NS-FOV(1)+s*NY,:);
         end
         xb(1:FOV(1)-oFRedT,:)=xb(1:FOV(1)-oFRedT,:)+x(NS-FOV(1)+1+oFRedT:NS,:);
         xb(NY-FOV(2)+1+oFRedT:NY,:)=xb(NY-FOV(2)+1+oFRedT:NY,:)+x(1:FOV(2)-oFRedT,:);
         xb=reshape(xb,[NY N(2:end)]);
         
         if m~=1
             xb=permute(xb,perm);
         end    
    %else  
    %     N=size(x);        
    %     x=reshape(x,[N(1:m-1) NY ceil(NS/NY) prod(N(m+1:nDimsIn))]);                    
    %     %if m==1%Bug in GPU when summing through the second dimension for
    %     an array of size 74x2x409600. Potentially related with hitting
    %     the memory limits
    %     %    x=permute(x,[1 3 2]);
    %     %    xb=sum(x,3);
    %     %else
    %     xb=sum(x,m+1);
    %     %end
    %     xb=reshape(xb,[N(1:m-1) NY N(m+1:nDimsIn)]);        
    %     perm=zeros(1,nDimsIn);perm(m)=-FOV(2);         
    %     xb=circshift(xb,perm);
    % end
else
    xb=x;
end
