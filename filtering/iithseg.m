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
info = dicominfo(fullfile(filefolder,filenames{1}))

voxel_size = [info.PixelSpacing; info.SliceThickness]';

I         = dicomread(fullfile(filefolder, filenames{1}));
classI    = class(I);
sizeI     = size(I);
numImages = length(filenames);
toc

%%

tic
hWaitBar = waitbar(0,'Reading dicom images');

ct = zeros(sizeI(1), sizeI(2), numImages, classI);

for i=100:-1:1
    fname = fullfile(filefolder,filenames{i});
    ct(:,:,i) = uint16(dicomread(fname));
    
    waitbar((length(filenames)-i+1)/length(filenames))
end

delete(hWaitBar)

volshow(ct)
toc

%%

tic
v = dicomreadVolume('/home/biom3d/akkm/IITH/S4_1001_1200/');
V = squeeze(v);

V(V<10) = 0;
% V(V>11) = 0;
imbinarize(V);
volshow(V)
niftiwrite(V,'z2_chloro_1000-1200');

toc

%% %% NiFTI file reduction
clear;
clc;

%ni2 = niftiread("D:\akkm\matlab\seg\reduction\z2_chloro.nii");
ni2 = niftiread("D:\akkm\matlab\seg\reduction\palak_z3.nii");
[x,y,z] = size(ni2);
%writematrix(ni2, 'chloro2_orig.txt');
writematrix(ni2, 'palak3_orig.txt');
num=0;
c=0;
sum_thresh = 0;

%% Extracting one slice
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

%%
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

%%
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

%%
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

%%
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

%%
ni2_new = ni2(temp_x0:temp_x1, temp_y0:temp_y1, :);

% volshow(ni2, 'Colormap',colormap)
% figure
% hold on
% volshow(ni2_new);

%writematrix(ni2_new, 'chloro2_reduced.txt');
writematrix(ni2_new, 'palak3_reduced.txt');
niftiwrite(ni2_new,'palak3_reduced');


