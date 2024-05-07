%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Sébastien Tosi (sebastien.tosi@gmail.com)
% Version: 1.1 (transwell, 0.5 XY downsampling from article conditions)
% Date:    06/03/24
%
% Input: Spinning disk 3D stack (3 channels): 
%        nuclei marker, primary cilia marker, cilia protein of interest 
%
% Assumptions:
%        - Primary cilium:
%              * overlaping with associated nucleus
%              * brightest object for cilia marker inside nucleus 
%              * volume within user bounds
%              * extending outward (base closest point to nucleus centroid)
%        - For accurate (and consistent) cilia length estimation, the image 
%          ZRatio should be close to 1 (Z step = pixel size). 
%          If this is critical to your study, resample the images
%          accordingly.
%
% Usage:
%       1) Set the image name (variable ImgFile)
%          If the image is not in the script folder, write its full path 
%       2) Set the imaging settings according to the acquisition conditions
%       3) Set the main measurement from section --> Main measurement <--
%          The default is the average intensity of the cilium protein
%       4) Run the script (green arrow)
%       5) If the nuclei/cilia detection are inaccurate, adjust sample 
%          parameters, especially: NucThr, CilThr, MinCilInt. Run again.
%          To make it easier, you may initially set MinCilVol = 0
%          and MaxCilVol = Inf to keep all detected primay cilia.
%       6) Fine tune the additional parameters if needed. These are
%          geometry only dependent and should not vary much between
%          samples (here adjusted for the pixel size reported in the 
%          article and after downsampling the images by 0.5 in XY).
%
% Output: 
%        Images (maximum intensity projections):
%        1) Segmented nuclei (if ShowNucImg == 1)
%        2) Cilia candidates (thin red), primary cilia (green),
%           Nuclei with (magenta) and without (red) primary cilium,
%           Nuclei centroids (blue), primary cilia bases (yellow)
%        3) Nuclei with valid primary cilium (blue),
%           color-coded main measurement centered on cilia bases
%
%        Matlab log:
%        - Maximum intensity in nuclei and primary cilia marker channels
%        - Total number of nuclei and nuclei with valid primary cilium
%        - Primary cilium base, body and average intensity (cilia protein)
%        - Primary cilium length (pixels) and base to body intensity ratio
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
iptsetpref('ImshowBorder','tight');
clear all;
close all;
clc;

%% Image File (XY downsampled by 0.5 for 2048x2048 image format)
ImgFile = 'Image_10371_0_5_crop.tif';

%% Imaging Settings
NChan = 3;              % Number of channels
CilChan = 1;            % Cilia marker channel index    
MeasChan = 2;           % Cilia measurement channel
NucChan = 3;            % Nuclei channel index
ZRatio = 2;             % Image ZRatio (Zstep / pixel size)
%% Sample Parameters
NucThr = 42;            % Nuclei detection intensity threshold (nuclei marker)
CilThr = 50;            % Cilia detection intensity threshold (cilia marker)
MinCilInt = 50;         % Primary cilia: minimum average intensity (cilia marker) 
MinCilVol = 50;         % Primary cilia: minimum volume (voxels) 
MaxCilVol = 875;        % Primary cilia: maximum volume (voxels)
%% Additional Parameters
NucGaussRad = 10;       % Nuclei channel: Gaussian radius (pixels)
NucBckRad = 70;         % Nuclei channel: background subtraction radius (pixels)
NucSplitRad = 12;       % Nuclei channel: splitting radius (pixels)
CilGaussRad = 1;        % Cilia channel: Gaussian radius (pixels) 
CilBox = [18 18 12];    % Cilia channel: local theshold box (pixels)
CilBaseMeasDst = 2;     % Primary cilia: base extension (voxels)
BckCorrBox = [128 128 1]; % Protein chan: background correction (pixels)
%% Display only
ShowNucImg = 1;         % Show nuclei channel + segmentation overlay
DilCil = 1;             % Cilium dilation (pixels)
ShowPrimCil = 1;        % Show primary cilia overlay (green)
SatCil = 1000;          % Cilium marker: intensity saturation level
SatMeas = 250;          % Cilium protein: intensity saturation level

%% Retrieve image Information
info = imfinfo(ImgFile);
NZ = numel(info)/NChan;
W = info(1).Width;
H = info(1).Height;
ZRef = round(NZ/2);

%% Load Image Stack
A = zeros(H,W,NZ,NChan,'uint16');
for z = 1:NZ
    for c = 1:NChan
        A(:,:,z,c) = imread(ImgFile,c+(z-1)*NChan);
    end
end

%% Force odd size box
CilBox = 2*floor(CilBox/2)+1;
if(BckCorrBox > 0)
    BckCorrBoxDs = 2*floor(0.25*BckCorrBox/2)+1;
else 
    BckCorrBoxDs = 0;
end

%% Nuclei 2D Segmentation (from nuclei channel Z-projection)
NucImg = A(:,:,:,NucChan);
PNucImg = max(NucImg,[],3);
disp(['Nuclei channel maximum intensity  : ' num2str(max(PNucImg(:)))]);
% Z-project nuclei channel + prefilter + threshold
PNucFlt = imgaussfilt(imtophat(PNucImg,strel('disk',NucBckRad)),NucGaussRad);
PNucMsk =  PNucFlt >= NucThr;
% Split touching nuclei (binary watershed) and remove nuclei touching image borders
PNucAll = PNucMsk & watershed(imgaussfilt(-bwdist(~PNucMsk,'euclidean'),NucSplitRad));
PNucMsk = imclearborder(PNucAll);
% Analyze connected components + build associated label mask
NucPixInds = bwconncomp(PNucMsk).PixelIdxList;
PNucLbl = zeros(size(PNucMsk));
for i = 1:numel(NucPixInds)
    PNucLbl(NucPixInds{i}) = i;
end
% Find nuclei centroids
NucCC = reshape([regionprops(PNucMsk,'Centroid').Centroid],2,[]).';
NNuc = size(NucCC,1);
% Instantiate result table
NucPrimCil = zeros(NNuc,9);

%% Cilia 3D Segmentation (from cilia marker channel)
CilImg = A(:,:,:,CilChan);
disp(['Pimary cilia channel maximum intensity : ' num2str(max(CilImg(:)))]);
% Pre-filter + local mean thresholding
CilCandMsk = uint8(imgaussfilt(CilImg,CilGaussRad)>(imboxfilt3(CilImg,CilBox)+CilThr));
% Analyze connected components
ObjsPixInds = bwconncomp(CilCandMsk).PixelIdxList;

%% Identify nuclei with valid primary cilium (min/max volume + min avg. 
%%   intensity + brightest object majoratirly overlaping nucleus)
% Loop over all cilia and fill first two columns of NucPrimCil: 
% Primary cilium label + average intensity
for i = 1:numel(ObjsPixInds)
    ObjPixInds = ObjsPixInds{i};
    [Y, X, Z] = ind2sub(size(CilImg), ObjPixInds);
    % Identify associated nucleus label (largest overlap)
    Inds = sub2ind(size(PNucMsk),Y,X);
    NucCands = PNucLbl(Inds);
    NucLbl = mode(NucCands(NucCands>0));
    % If cilium majoritarily overlaps with a nucleus 
    if NucLbl > 0
        ObjAvgInt = mean(CilImg(ObjPixInds));
        % If cilium checks all conditions
        if ObjAvgInt > MinCilInt && numel(ObjPixInds) >= MinCilVol && numel(ObjPixInds) < MaxCilVol
            % If cilium is brigher than previous valid primary cilium
            if ObjAvgInt > NucPrimCil(NucLbl,2)
                NucPrimCil(NucLbl, 1) = i;
                NucPrimCil(NucLbl, 2) = ObjAvgInt;
            end
        end
    end
end

%% Identify primary cilia base locations (closest cilium point to nucleus centroid)
%% report and create a mask holding primary cilia + bases + nuclei with primary cilium
% Fill next three columns of NucPrimCil: Base CX, Base CY, Base CZ
CilMsk = uint8(zeros(size(CilCandMsk)));
PNucCilMsk = uint8(zeros(size(PNucMsk)));
CilBaseInds = [];
for i = 1:NNuc
    IndCilCand = round(NucPrimCil(i,1));
    if IndCilCand > 0
        PNucCilMsk(NucPixInds{i}) = 255;            % Highlight associated nucleus
        CilMsk(ObjsPixInds{IndCilCand}) = 255;      % Highlight primary cilium
        CilMarkerAvgInt = mean(CilImg(ObjsPixInds{IndCilCand}));
        % Find base (minimum distance to nucleus centroid)
        [Y, X, Z] = ind2sub(size(CilImg), ObjsPixInds{IndCilCand});
        [val, ind] = min((X-NucCC(i,1)).^2 + (Y-NucCC(i,2)).^2 + ((Z-ZRef)*ZRatio).^2);
        NucPrimCil(i, 3:5) = [X(ind) Y(ind) Z(ind)];
        CilBaseInds = [CilBaseInds sub2ind(size(CilMsk),Y(ind),X(ind),Z(ind))];
    end
end

%% Find the indices of nuclei with a primary cilium
inds = find(any(NucPrimCil(:,3:5),2));
NPrimCil = numel(inds);
disp(['Total number of nuclei     : ' num2str(NNuc)]);
disp(['Nuclei with primary cilium : ' num2str(NPrimCil)]);

%% Geodesic tracing from cilia base start to base extension
D = bwdistgeodesic(logical(CilCandMsk), CilBaseInds);

%% Measure signal inside cilia (average and base)
% NucPrimCil last 4 columns:
% Cil Base avg. MeasInt, Cil Body avg. MeasInt, Cil avg. MeasInt, Cil Lgth.
MeasImg = A(:,:,:,MeasChan);
% Correct local background
if sum(BckCorrBoxDs)>0
    BckEst = imresize3(medfilt3(imresize3(MeasImg, 0.25), BckCorrBoxDs), size(MeasImg));
    MeasImg = MeasImg - BckEst;
end
PrimCilBaseCoords = NucPrimCil(:, 3:5);
for i = 1:NPrimCil
    CilPixInds = ObjsPixInds{NucPrimCil(inds(i),1)};
    CilBasePixInds = CilPixInds(D(CilPixInds)<=CilBaseMeasDst);
    CilBodyPixInds = CilPixInds(D(CilPixInds)>CilBaseMeasDst);
    NucPrimCil(inds(i),6) = mean(MeasImg(CilBasePixInds));
    NucPrimCil(inds(i),7) = mean(MeasImg(CilBodyPixInds));
    NucPrimCil(inds(i),8) = mean(MeasImg(CilPixInds));
    NucPrimCil(inds(i),9) = max(D(CilPixInds));
    % For short cilia, set ratio to 1 (base / body split is arbitrary)
    if(NucPrimCil(inds(i),6) == 0 || NucPrimCil(inds(i),7) == 0)
        NucPrimCil(inds(i),6) = 1e-6;
        NucPrimCil(inds(i),7) = 1e-6;
    end
end
disp('   Base Int  Body Int  Avg Int   Length     Base/Body Int');
Ratios = NucPrimCil(:,6)./NucPrimCil(:,7);
disp([NucPrimCil(inds,6) NucPrimCil(inds,7) NucPrimCil(inds,8) NucPrimCil(inds,9) Ratios(inds)]);

%% Create result images
PNucFlt8b = uint8(single(PNucFlt)/single(max(PNucFlt(:)))*255);
PNucAllEdg = 255*uint8(edge(PNucAll));
PNucEdg = 255*uint8(edge(PNucCilMsk));
PCil8b = uint8((255*single(max(CilImg,[],3))/SatCil));
PCilCandEdg = 255*uint8(edge(imdilate(max(CilCandMsk,[],3),strel('disk',1))));
PCilEdg = 255*uint8(edge(imdilate(max(CilMsk,[],3),strel('disk',DilCil))));
CilAvgMeasInt = NucPrimCil(:,8);
PNucAllEdg = imdilate(PNucAllEdg, strel(ones(3,3)));
PNucEdg = imdilate(PNucEdg, strel(ones(3,3)));

%% Display nuclei channel (optional)
if ShowNucImg
    Mrg = uint8(cat(3,PNucFlt8b,PNucFlt8b,PNucFlt8b+PNucAllEdg));
    figure('MenuBar','none','toolbar','figure', 'Name', 'Nuclei');
    imshow(Mrg);hold on;axis equal;axis off;
end

%% Display cilia channel with nuclei centroids and cilia bases
Mrg = uint8(cat(3,PCil8b+PNucAllEdg+PCilCandEdg/4,PCil8b+ShowPrimCil*PCilEdg,PCil8b+PNucEdg));
figure('MenuBar','none','toolbar','figure', 'Name', 'Cilia marker (green: primary cilia with yellow base)');
imshow(Mrg);hold on;axis equal;axis off;
plot(NucCC(inds,1),NucCC(inds,2),'bx');axis off;
plot(NucPrimCil(inds,3),NucPrimCil(inds,4),'yx');

%% --> Main measurement <--
% From table NucPrimCil primary cilia measurements (columns 6 to 9): 
% 6: Base meas. average intensity, 7: Body meas. average intensity
% 8: Cilim average intensity, 9: Cilium length (pixels)
Meas = NucPrimCil(:,8);
% To estimate the ratio between the intensity of the cilia base and body:
%Meas = NucPrimCil(:,6)./NucPrimCil(:,7);

%% Display protein channel and main measurement (color-coded)
ProjMeas = single(max(MeasImg,[],3));
ProjMeas8b = uint8((ProjMeas/SatMeas)*255);
Mrg = uint8(cat(3,ProjMeas8b,ProjMeas8b+0*PCilEdg,ProjMeas8b+PNucEdg));
Cols = hot(255);
figure('MenuBar','none','toolbar','figure', 'Name', 'Cilia protein (color-coded main measurement)');
imshow(Mrg, [0 255]);colormap('gray');hold on;
% Normalize color-coding to 8-bit range
NmrMeas = 55+round(200*(Meas-min(Meas(inds)))./(max(Meas(inds))-min(Meas(inds))));
scatter(NucPrimCil(inds,3),NucPrimCil(inds,4),250,Cols(NmrMeas(inds),:),'o');
scatter(NucPrimCil(inds,3),NucPrimCil(inds,4),25,Cols(NmrMeas(inds),:),'x');
