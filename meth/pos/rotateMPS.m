function x=rotateMPS(x,ARot)

%ROTATEMPS rotates the MPS frame
%   X=ROTATEMPS(X,AROT) applies a rotation from the MPS reference frame to
%   the ijk reference frame
%   X is the image to be rotated
%   AROT is a rotation matrix
%   It returns:
%   X, the rotated image
%

for m=1:3
    s(m)=find(ARot(m,:));
end
s(4)=4;
x=permute(x,s);
for m=1:3
    v(m)=ARot(m,s(m));
end
if v(1)<0;x=x(end:-1:1,:,:,:);end
if v(2)<0;x=x(:,end:-1:1,:,:);end
if v(3)<0;x=x(:,:,end:-1:1,:);end
