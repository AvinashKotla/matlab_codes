clear;
clc;

%%

tic
filefolder = fullfile('/home/biom3d/akkm/IITH/S4dcm/');
files = dir(fullfile(filefolder, '*.dcm'));
filenames = {files.name};

info = dicominfo(fullfile(filefolder,filenames{1}));
I         = dicomread(fullfile(filefolder, filenames{1}));
classI    = class(I);
sizeI     = size(I);
numImages = length(filenames);

zone_size = 300;
zones = ceil(numImages/zone_size);

v = dicomreadVolume('/home/biom3d/akkm/IITH/S4dcm/');
ct = squeeze(v);

toc

%%

tic

% ct1 = ct(:,:,zone_size*0+1:zone_size*1); % lower=12, upper=25
% ct2 = ct(:,:,zone_size*1+1:zone_size*2);
% ct3 = ct(:,:,zone_size*2+1:zone_size*3);
% ct4 = ct(:,:,zone_size*3+1:zone_size*4);
% ct5 = ct(:,:,zone_size*4+1:zone_size*5);
ct6 = ct(:,:,zone_size*5+1:zone_size*6);
% ct7 = ct(:,:,zone_size*6+1:zone_size*7);
% ct8 = ct(:,:,zone_size*7+1:zone_size*8);
% ct9 = ct(:,:,zone_size*8+1:zone_size*9);
% ct10 = ct(:,:,zone_size*9+1:zone_size*10);
% ct11 = ct(:,:,zone_size*10+1:zone_size*11);
% ct12 = ct(:,:,zone_size*11+1:end);

toc

%% %% NiFTI file reduction

tic

ni2 = ct(:,:,2500:2600);
ni22 = ni2;

[x,y,z] = size(ni2);
num=0;
c=0;
sum_thresh = 1;

ni2(ni2 < 13) = 0;
ni2(ni2 > 25) = 0;
% volshow(ni2)
toc

%%
tic
aa = ct(:,:,101:200);
aa(aa<13) = 0;
aa(aa>25) = 0;

toc

%% Extracting one slice

tic
ni2d = ni2(:,:,26);
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
% to compute the upper bounding box horizontal line
tic
temp_x0 = x;
for k=1:z % for all slices
    for i=1:x % starting from top
%         fprintf("%d\n",i)
        if(sum(ni2(i,:,k)) > sum_thresh) % not less than, maximum of non zero sum rows
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
% to compute the lower bounding box horizontal line
tic
temp_x1 = 0;
arr = zeros(1,z);

for k=1:z
    for i=x:-1:1 % starting from bottom
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
% to compute the left bounding box vertical line
tic
temp_y0 = y;
for k=1:z
    for i=1:y % starting from left
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
% to compute the right bounding box vertical line
tic
temp_y1 = 0;
for k=1:z
    for i=y:-1:1 % starting from right
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

volshow(ni22, 'Colormap',colormap)
figure
hold on
volshow(ni2_new,'Colormap',colormap);

toc

%%

tic

save ni2new.mat ni2_new;

toc

%%

tic

load ni2new.mat;


toc