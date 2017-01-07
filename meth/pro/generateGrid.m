function rGrid=generateGrid(N,fo,Nres,cent,nor,gpu,rGrid0)

%GENERATEGRID   Generates spatial/k-spatial grids
%   KGRID=GENERATEGRID(N,FO,NRES,CENT,NOR,GPU,RGRID0) generates a grid for 
%   a space of size N
%   N are the dimensions of the space
%   FO determines whether the grid corresponds to k-space (1) or space
%   (0, default)
%   NRES defines the size of the output grid
%   CENT defines the center of coordinates (empty, default)
%   NOR indicates wheter to normalize the space (default 1)
%   GPU is a flag that determines whether to generate gpu (1) or cpu (0) 
%   arrays
%   rGrid0 is a precomputed grid, defaults to empty
%   It returns:
%   KGRID, the generated k-spatial grid, as a cell with each entry 
%   specifying a given direction of the space
%

if ~exist('fo','var');fo=0;end
if ~exist('Nres','var');Nres=N;end
if ~exist('cent','var');cent=[];end
if ~exist('nor','var');nor=1;end
if ~exist('gpu','var');gpu=0;end

NDims=length(N);
res=N./Nres;
rGrid=cell(1,NDims);
for m=1:NDims    
    if isempty(cent)
        if fo==0
            if ~exist('rGrid0','var')
                rGrid{m}=-floor(N(m)/2)+0.5*(res(m)-1):res(m):ceil(N(m)/2)-1;
            else
                rGrid{m}=rGrid0{m}(1)+0.5*(res(m)-1):res(m):rGrid0{m}(end);
            end
        else
            rGrid{m}=-floor(Nres(m)/2):ceil(Nres(m)/2)-1;
        end
    else
        if fo==0
            rGrid{m}=(1+0.5*(res(m)-1):res(m):N(m))-cent(m);
        else
            rGrid{m}=(1:Nres(m))-cent(m);
        end
    end
    perm=1:NDims;perm(2)=m;perm(m)=2;
    rGrid{m}=single(permute(rGrid{m},perm));
    if gpu>0
        rGrid{m}=gpuArray(rGrid{m});
    end
    if nor || fo~=0
        rGrid{m}=rGrid{m}/N(m);
    end
    if fo~=0
        rGrid{m}=fo*2*pi*rGrid{m};
    end
end
