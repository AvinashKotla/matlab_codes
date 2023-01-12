%clc;
clear;
clc;

%ni3 = niftiread("F:\Work\MATLAB\seg\z2_chloro.nii");
%[x,y,z] = size(ni3);

%%
filefolder = fullfile('/home/biom3d/akkm/IITH/S4/Recon/');
% filefolder = fullfile('F:\Work\IISER_Pune_CT\chlorophyta\chloro_highres\chlorophyta_19_20220514_145135');
% filefolder = fullfile('F:\Work\IISER_Pune_CT\series-000001');
files = dir(fullfile(filefolder, '*.dcm'));
filenames = {files.name};


%%

info = dicominfo(fullfile(filefolder,filenames{1}))

voxel_size = [info.PixelSpacing; info.SliceThickness]';

I         = dicomread(fullfile(filefolder, filenames{1}));
classI    = class(I);
sizeI     = size(I);
numImages = length(filenames);

%%
hWaitBar = waitbar(0,'Reading dicom images');

ct = zeros(sizeI(1), sizeI(2), numImages, classI);

for i=length(filenames):-1:1
    fname = fullfile(filefolder,filenames{i});
    ct(:,:,i) = uint16(dicomread(fname));
    
    waitbar((length(filenames)-i+1)/length(filenames))
end

delete(hWaitBar)

%%
ct = flipdim(ct,3);
im = ct(:,:,250);
max_level = double(max(im(:)));
imt = imtool(im, [0,max_level]);
% montage({im, imadjust(ct(:,:,250))});

%%
z = 150;
ctTemp = ct(:,:,z);
% imto = imtool(ctTemp);
lb = 3200;
ub = 4000;
% lb = 850;
% ub = 1200;


ctTemp(ctTemp <= lb) = 0;
ctTemp(ctTemp >= ub) = 0;
% ctTemp(end-80:end,:) = 0;

imA = imadjust(ctTemp);
imshow(imA);

CC = bwconncomp(im,4)
reg = regionprops(CC,'all')


%%
slice = ct(:,:,200);
slice2 = ct(:,:,260);
imt = imtool(slice);

%%
BW = imcomplement(BW);
BW = imclearborder(BW);
BW = imfill(BW, "holes");
radius = 3;
decomposition = 0;
se = strel("disk",radius,decomposition);
BW = imerode(BW, se);
maskedImageXY = XY;
maskedImageXY(~BW) = 0;
imshow(maskedImageXY)



%%
ctadjust = ct;
ctadjust(ctadjust <= lb) = 0;
ctadjust(ctadjust >= ub) = 0;
ctadjust(end-80:end,:,:) = 0;

bw = ctadjust > 0;
nhood = ones([7 7 3]);

bw = imopen(bw,nhood);
imshow(bw(:,:,z));

%%
L = bwlabeln(bw);
stats = regionprops(L,'Area','Centroid');

LL = L(:,:,z) + 1;
cmap = hsv(length(stats)); 
cmap = [0 0 0; cmap];
LL = cmap(LL,:);
LL = reshape(LL,[sizeI,3]);

imshow(LL);

%%
A = [stats.Area];
biggest = find(A == max(A));
ctadjust(L ~= biggest) = 0;

imA = imadjust(ctadjust(:,:,z));
imshow(imA);

%%
level = thresh_tool(unint16(ctadjust(:,:,z)),'gray');

ctpartition = uint8(zeros(size(ctadjust)));
ctpartition(ctadjust < level & ctadjust > 0) = 2;
ctpartition(ctadjust >= level) = 3;

imshow(ctpartition(:,:,z),[0 0 0;0 0 0;0.25 0.25 1;1 1 1])

%% 
% tmp = text(128,200,'Rendering: Wait','color','r','fontsize',1

