%MULTISLICEALIGNEDSENSE_EXP1 script performs the experiment included in
%Figure 3 of the manuscript ''3D motion corrected SENSE reconstruction for
%multi-shot multi-slice MR: application to neonatal brain imaging'', L.
%Cordero-Grande, E. J. Hughes, J. Hutter, A. N. Price, and J. V. Hajnal
addpath(genpath('.'));
pathOu='./data';%Data path
gpu=1;%0->Use CPU / 1->Use GPU
debug=1;%0/1/2 increasing the amount of debug info provided 

%Load synthetic data:
% - Measured data y of size PE-readout-slice-coil
% - Sensitivity maps S of size PE-readout-slice-coil
% - Reconstruction mask W of size PE-readout-slice
% - Sampling trajectory Ak of size 
% - Slice thickness SlTh
% - Slice overlap SlOv
% - Encoding values Encoding
% - Rotation matrix Rot
% - Image matrix Phi

%Code for uncompressed data
%nameX=sprintf('%s/yT2',pathOu);load(nameX);
%perc=0.95;
%[S,y]=coilArrayCompression(S,y,perc,gpu);
%save(sprintf('%s/yT2Comp',pathOu),'y','S','W','Ak','SlTh','SlOv','Encoding','Rot','Phi');
%return

%Load compressed data
nameX=sprintf('%s/yT2Comp',pathOu);load(nameX);

%Parameters fixed in this high level implementation
parX.toler=1e-6;parX.gibbsRing=[0.4 0.4];parX.alpha=[40 20];%For T2MS


Vpar.threeD=[0 1 1 1 1];
Vpar.outl=[0 0 0 1 1];
Vpar.correct=[0 0 1 1 1];
Vpar.thplc=[0 0 1 0 1];

%Initialization
NA=size(Ak);NX=size(W);
T=single(zeros([1 1 1 1 NA(5) NX(3) 6]));
outlD=single(ones([1 1 NX(3) 1 NA(5)]));
x=zeros(NX);
mNorm=max(abs(y(:)));
NS=size(S);NY=size(y);

%recType=1: NMC
%recType=2: NMC-SP 
%recType=3: MC-NOU
%recType=4: MC-NTP
%recType=5: MC
suffix={'NMC','NMC-SP','MC-NOU','MC-NTP','MC'};
suffixFull={'Conventional uncorrected SENSE (NMC)','Uncorrected with slice profile filter (NMC-SP)','Corrected without outlier rejection (MC-NOU)','Corrected without through-plane motion (MC-NTP)','Fully corrected (MC)'};
for recType=1:length(Vpar.threeD)
    if debug>0;fprintf('\nReconstructing %s...\n',suffixFull{recType});trec=tic;end
    %Initialization    
    x(:)=0;T(:)=0;outlD(:)=1;
    parX.threeD=Vpar.threeD(recType);parX.outl=Vpar.outl(recType);parX.correct=Vpar.correct(recType);parX.thplc=Vpar.thplc(recType);
    %2 multirresolution levels L in case of motion correction
    if parX.correct==0;estT=0;else estT=[1 0];end
    if ~parX.outl;parX.outlP=inf;else parX.outlP=1.2;end
    L=length(estT);
    H=sliceProfile(NX(3),SlTh,SlOv,parX);

    for l=1:L
        if debug==2;fprintf('Resolution level: %d\n',l);tstart=tic;end
        res=2^(L-l);
        NXres(3)=NX(3);NYres(3)=NY(3);
        NXres(1:2)=floor(NX(1:2)/res);NYres(1:2)=floor(NY(1:2)/res);        
        xRes=resampling(x,NXres);SRes=resampling(S,NXres);
        yRes=resampling(y,NYres);
        WRes=resampling(W,NXres);
        WRes=single(abs(WRes)>0.5);    
        AkRes=resampling(Ak,NYres(1),1);
        cent=[1 1 1];
        xGrid=generateGrid(NX,0,NXres,cent,0); 
        if debug==2;telapsed=toc(tstart);fprintf('Time setting multires: %.4fs\n',telapsed);end
        
        [xRes,T,outlD]=optimizeLevelMS2D(xRes,yRes,T,SRes,WRes,AkRes,xGrid,mNorm,parX,debug,res,estT(l),gpu,SlTh,SlOv,outlD,parX.alpha(l));  
        x=W.*resampling(xRes,NX); 
    end
    xF=permute(x,[2 1 3]);WF=permute(W,[2 1 3]);
    %To perform superresolution something along this lines:
    %NDims=3;
    %parX.gibbsRing(3)=parX.gibbsRing(2); 
    xF=filtering(xF,H);    
    NDims=2;
    
    %Gibbs ringing
    xF=gibbsRingingFilter(xF,NDims,parX.gibbsRing);

    %Zero filling
    N=[];
    N(1)=Encoding.XRes;
    N(2)=round(Encoding.YRes*Encoding.KyOversampling);
    WF=resampling(WF,N);
    xF=single(abs(WF)>0.5).*resampling(xF,N);

    %Remove oversampling
    pady=(size(xF,2)-Encoding.YRes)/2;
    padya=floor(pady);padyb=ceil(pady);
    xF=xF(:,1+padya:end-padyb,:);

    %Zero fill
    padi=[Encoding.XReconRes-Encoding.XRes Encoding.YReconRes-Encoding.YRes]/2;
    padia=floor(padi);padib=ceil(padi);
    xF=padarray(xF,padia,0,'pre');xF=padarray(xF,padib,0,'post');

    %Rotate image
    Rot=Rot(1:3,1:3);
    xF=rotateMPS(xF,Rot);
 
    %Write to nifti file
    writeNIFTI(abs(xF),Phi,pathOu,'xT2',suffix{recType});
    if debug>0;telapsed=toc(trec);fprintf('Time reconstructing %s: %.2fs\n',suffixFull{recType},telapsed);end
end

