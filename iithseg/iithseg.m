clear;
clc;


%%

tic
filefolder = fullfile('/home/biom3d/akkm/IITH/S3/Recon_2/');
% filefolder = fullfile('F:\Work\IISER_Pune_CT\chlorophyta\chloro_highres\chlorophyta_19_20220514_145135');
% filefolder = fullfile('F:\Work\IISER_Pune_CT\series-000001');
files = dir(fullfile(filefolder, '*.dcm'));
filenames = {files.name};
toc

%%

tic
info = dicominfo(fullfile(filefolder,filenames{1}));

voxel_size = [info.PixelSpacing; info.SliceThickness]';

I         = dicomread(fullfile(filefolder, filenames{1}));
classI    = class(I);
sizeI     = size(I);
numImages = length(filenames);
toc

%%

zone_size = 300;
zones = ceil(numImages/zone_size);
% s.ct = zeros(1,zones);
% s.ct = zeros(sizeI(1), sizeI(2), numImages, classI);

%%

tic
hWaitBar = waitbar(0,'Reading dicom images');

ct = zeros(sizeI(1), sizeI(2), numImages, classI);
% ctzone = zeros(sizeI(1), sizeI(2), numImages, classI);
i=0;
j=1;
% for j=1:zones
for i=zone_size*i+1:1:zone_size*j
    fname = fullfile(filefolder,filenames{i});
    ct(:,:,i) = uint16(dicomread(fname));
    waitbar((length(filenames)-i+1)/length(filenames))
end
%     j=j+1;
% end

delete(hWaitBar)

% volshow(ct)
toc

%%
tic

ct1 = ct(:,:,zone_size*0+1:zone_size*1);
ct2 = ct(:,:,zone_size*1+1:zone_size*2);
ct3 = ct(:,:,zone_size*2+1:zone_size*3);
ct4 = ct(:,:,zone_size*3+1:zone_size*4);
ct5 = ct(:,:,zone_size*4+1:zone_size*5);
ct6 = ct(:,:,zone_size*5+1:zone_size*6);
ct7 = ct(:,:,zone_size*6+1:zone_size*7);
ct8 = ct(:,:,zone_size*7+1:zone_size*8);
ct9 = ct(:,:,zone_size*8+1:zone_size*9);
ct10 = ct(:,:,zone_size*9+1:zone_size*10);
ct11 = ct(:,:,zone_size*10+1:zone_size*11);
ct12 = ct(:,:,zone_size*11+1:end);

toc

%%
tic

ct12(ct12<10) = 0;
ct12(ct12>11) = 0;


% volshow(ct1);


toc

%%

tic
% v = dicomreadVolume('/home/biom3d/akkm/IITH/S3dcm_501-800/');
v = dicomreadVolume('/home/biom3d/akkm/IITH/S3dcm/');
% v = dicomreadVolume('/home/biom3d/akkm/IITH/S3dcm_501-800/');
V = squeeze(v);

% V(V<10) = 0;
% V(V>11) = 0;
% imbinarize(V);
% volshow(V)

toc

%% %% NiFTI file reduction

tic
%ni2 = niftiread("D:\akkm\matlab\seg\reduction\z2_chloro.nii");
ni2 = niftiread("/home/biom3d/akkm/medical_imaging/filtering/iithseg/z2_chloro_1000-1200.nii");
[x,y,z] = size(ni2);
%writematrix(ni2, 'chloro2_orig.txt');
writematrix(ni2, 'z2_chloro_1000-1200.txt');
num=0;
c=0;
sum_thresh = 0;
toc

%% Extracting one slice
tic
ni2d = ni2(:,:,78);
[x2d,y2d] = size(ni2d);
fprintf("Before: %d, %d\n",size(ni2d));

for i=1:x
    for j=1:y
        if (ni2d(i,j) ~=0)
            fprintf("%d, %d, %d\n",i,j,c);
            return
        end
    end
end
toc

%%
tic
temp_x0 = x;
for k=1:z
    for i=1:x
%         fprintf("%d\n",i)
        if(sum(ni2(i,:,k)) > sum_thresh)
%             fprintf("%d, %d, %d\n",i,j,c);
            if i < temp_x0
                num = num+1;
                temp_x0 = i;
                fprintf('%d, %d\n',temp_x0,k)
                break
            end
        end
    end
end
toc
%%
tic
temp_x1 = 0;
arr = zeros(1,z);

for k=1:z
    for i=x:-1:1
%         fprintf("%d\n",i)
        if(sum(ni2(i,:,k)) > sum_thresh)
%             fprintf("%d, %d, %d\n",i,j,c);
            if i > temp_x1
                arr(num+1) = temp_x1;
                num = num+1;
                temp_x1 = i;
                fprintf('%d, %d\n',temp_x1,k)
                break
            end
        end
    end
end
toc

%%
tic
temp_y0 = y;
for k=1:z
    for i=1:y
%         fprintf("%d\n",i)
        if(sum(ni2(:,i,k)) > sum_thresh)
%             fprintf("%d, %d, %d\n",i,j,c);
            if i < temp_y0
                num = num+1;
                temp_y0 = i;
                fprintf('%d, %d\n',temp_y0,k)
                break
            end
        end
    end
end
toc

%%
tic
temp_y1 = 0;
for k=1:z
    for i=y:-1:1
%         fprintf("%d\n",i)
        if(sum(ni2(:,i,k)) > sum_thresh)
%             fprintf("%d, %d, %d\n",i,j,c);
            if i > temp_y1
                num = num+1;
                temp_y1 = i;
                fprintf('%d, %d\n',temp_y1,k)
                break
            end
        end
    end
end
toc

%%
tic
ni2_new = ni2(temp_x0:temp_x1, temp_y0:temp_y1, :);

% volshow(ni2, 'Colormap',colormap)
% figure
% hold on
% volshow(ni2_new);

%writematrix(ni2_new, 'chloro2_reduced.txt');
writematrix(ni2_new, 'z2_chloro_1000-1200_reduced.txt');
niftiwrite(ni2_new,'z2_chloro_1000-1200_reduced');
toc

%%

ni2bin = dec2bin(ni2_new);