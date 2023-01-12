clear;
clc;

%%
tic
v = dicomreadVolume('F:\Work\IISER_Pune_CT\chlorophyta\chloro_highres\chlorophyta_19_20220514_145135');
V = squeeze(v);

XY = V(:,:,198);
XZ = squeeze(V(370,:,:));
toc
%%
imageSegmenter(XY)
%%
BW = imcomplement(BW);
se = strel('disk',8,0);
BW = imopen(BW,se);
BW = imclearborder(BW);
maskedXY = XY;
maskedXY(~BW) = 0;
% imshow(imadjust(maskedXY))

%%
imageSegmenter(XZ)

%%
tic
BW1 = imcomplement(BW1);
se = strel('disk',10,0);
BW1 = imopen(BW1,se);
BW1 = imclearborder(BW1);
maskedXZ = XZ;
maskedXZ(~BW1) = 0;
toc

%%
mask = false(size(V));
mask(:,:,198) = maskedXY;
mask(360,:,:) = mask(360,:,:)|reshape(maskedXZ,[1,512,427]);

%%
tic
V = histeq(V);
new = activecontour(V,mask,300,"Chan-Vese");
volshow(new)
toc

%%
tic
segmentedImage = V.*single(new);
toc