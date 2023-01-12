%% %% Edge detection (filter_1)

I = rgb2gray(imread("coronal_0427.png"));
im2 = imadjust(im);

subplot(2, 4, 1),
imshow(I);
title("Gray Scale Image");
 
% Sobel Edge Detection
J = edge(I, 'Sobel');
subplot(2, 4, 2),
imshow(J);
title("Sobel");
 
% Prewitt Edge detection
K = edge(I, 'Prewitt');
subplot(2, 4, 3),
imshow(K);
title("Prewitt");
 
% Robert Edge Detection
L = edge(I, 'Roberts');
subplot(2, 4, 4),
imshow(L);
title("Robert");
 
% Log Edge Detection
M = edge(I, 'log');
subplot(2, 4, 5),
imshow(M);
title("Log");
 
% Zerocross Edge Detection
M = edge(I, 'zerocross');
subplot(2, 4, 6),
imshow(M);
title("Zerocross");
 
% Canny Edge Detection
N = edge(I, 'Canny');
subplot(2, 4, 7),
imshow(N);
title("Canny");

subplot(2,4,8),
imshow(imadjust(I));
title("Adjusted");

%% %% DICOM Segmentation (readdicom)
%clc;
clear;
clc;

%ni3 = niftiread("F:\Work\MATLAB\seg\z2_chloro.nii");
%[x,y,z] = size(ni3);

%%
filefolder = fullfile(pwd, 'chloro_z3');
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
%imt = imtool(im, [0,max_level]);
montage({im, imadjust(ct(:,:,250))});

%%
imtool close all;

minct = min(ct(:));
maxct = max(ct(:));
montage(reshape(uint16(ct),[size(ct,1),size(ct,2),1,size(ct,3)]));
set(gca,'clim',[0,100]);

%%
z = 300;
ctTemp = ct(:,:,z);

lb = 2000;
ub = 5000;

ctTemp(ctTemp <= lb) = 0;
ctTemp(ctTemp >= ub) = 0;
ctTemp(end-80:end,:) = 0;

imA = imadjust(ctTemp);
imshow(imA);

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

%% Binary 3D model (binarymodel)

% kernel = zeros(3,3);
% k = [1 0 -1; 0 2 3; -1 3 -2];
% 
% new = conv2(im2, k);
% montage({im,new});

vol=zeros(512,512,100);
vol2 = vol;
slice = vol(:,:,1);

x0 = 255;
y0 = 255;
for z=1:100
    x0 = x0 + randi([-5 5],1);
    y0 = y0 + randi([-5 5],1);
    for x=1:512
        for y=1:512
            if 1.5*(x-x0)^2 + 1.5*(y-y0)^2  < 400 
                vol(x,y,z) = 255;
            end
        end
    end
end

% slice((x-x0)^2 + (y-y0)^2  < 100) = 255;
%imshow(slice);
volshow(vol);

%% %% Laplacian

x = 1:10;
% y = [1 3 4 6 -1 0 3 8 2 1];
y = randi([-10 10],[1 10]);
plot(x,y);

ynew = y;

for k=1:3
    for i=2:9
        ynew(i) = (ynew(i) + ynew(i+1))/2;
    end
    hold on
    plot(x,ynew);
end

%%
% x = [1 8];
% y = [2 10];
% 
% p = [1 8; 2 10];
% 
% plot(p(1,:),p(2,:));
% 
% hold on;
% theta = 30;
% ker = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% 
% pnew = ker * p;
% 
% plot(pnew(:,1),pnew(:,2));

%%

x = 1:10;
y = [1 3 4 6 -1 0 3 8 2 1];

p(:,1) = x;
p(:,2) = y;

xwidth = 300;
ywidth = 200;

figure (1)
h = figure(1);
set(gcf,'PaperPositionMode','auto')
set(h, 'Position', [0 0 xwidth ywidth])
plot(x,y)
saveas(gcf,'fig1.png')

% hold on;
% 
% % theta = 45 * pi/180;
% 
% theta = 0:9;
% i = theta * (pi/180);
% for i = 0:9
%     ker = [cos(i) -sin(i); sin(i) cos(i)] ;
%     % pnew = ker * p';
%     p(:,1) = x;
%     p(:,2) = y;
%     pnew = p * ker;
% 
%     % plot(pnew(1,:),pnew(2,:));
%     plot(pnew(:,1),pnew(:,2));
% end


%%
 % create the video writer with 1 fps
%  writerObj = VideoWriter('myVideo.avi');
%  writerObj.FrameRate = 1;
% 
%  % set the seconds per image
%  secsPerImage = [5 10 15];
% 
%  % open the video writer
%  open(writerObj);
% 
%  % write the frames to the video
%  for u=1:length(images)
%      % convert the image to a frame
%      frame = im2frame(images{u});
% 
%      for v=1:secsPerImage(u) 
%          writeVideo(writerObj, frame);
%      end
%  end
% 
%  % close the writer object
%  close(writerObj);

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

%% Branching
function [x,y] = branch(x0,y0,theta1, theta2,len)
% function branch(x0,y0,theta1, theta2,len)


xorig = x0 * ones(1,len+1);
yorig = [y0:1:y0+len];

x1 = xorig;
y1 = yorig;
x2 = xorig;
y2 = yorig;

% theta
theta1 = theta1*pi/180;
theta2 = theta2*pi/180;
rot1 = [cos(theta1) -sin(theta1); sin(theta1) cos(theta1)];
rot2 = [cos(-theta2) -sin(-theta2); sin(-theta2) cos(-theta2)];

new1 = [x1;y1];
new1 = rot1*new1;
x1 = new1(1,:);
y1 = new1(2,:);

new2 = [x2;y2];
new2 = rot2*new2;
x2 = new2(1,:);
y2 = new2(2,:);

x = [xorig x1 x2];
y = [yorig y1 y2];

% plot(xorig,yorig,x1,y1,x2,y2)
plot(x,y)

end


%% Fractal
function [x,y] = fractal(x1,y1,orientation,levels)

% x2 = x1+cosd(theta)*iterations;
% y2 = y1+sind(theta)*iterations;

x2 = x1+cosd(orientation)*levels;
y2 = y1+sind(orientation)*levels;
theta = 30;
if levels~=0
    line([x1 x2],[y1 y2],'LineWidth',2);
    x = x1:0.1:x2;
    y = y1:0.1:y2;
    fractal(x2,y2,orientation+theta,levels-1);
    fractal(x2,y2,orientation-theta,levels-1);
end

end


