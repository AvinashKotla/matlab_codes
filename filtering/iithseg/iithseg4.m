clear;
clc;

%% loading all the slices

tic
v = dicomreadVolume('/home/biom3d/akkm/IITH/S4dcm/');
ct = squeeze(v);
[height,width,numSlices] = size(ct);
zone_size = 100;
num_zones = ceil(numSlices/zone_size);
clear v;
toc

%%

tic
ct1 = ct(:,:,1501:1595);

ct1(ct1<12)=0;
ct1(ct1>25)=0;

toc;



%% file cropping

tic
clear v;
% x0 = zeros(1,height);
% x1 = zeros(1,height);
% y0 = zeros(1,height);
% y1 = zeros(1,height);
surfmax = zeros(1,num_zones-1);
ctsub = zeros(height,width,zone_size+1);
for q=0:5
   % tic
    if (zone_size*q+zone_size > numSlices)
        ctsub = ct(:,:,zone_size*q+1:end);
    else
        ctsub = ct(:,:,zone_size*q+1:zone_size*q+zone_size);
    end
    [x,y,z] = size(ctsub);
    
    num=0;
    c=0;
    
    lower_thresh = 13;
    upper_thresh = 25;
    sum_thresh = 1;
    

    ctsub(ctsub < lower_thresh) = 0;
    ctsub(ctsub > upper_thresh) = 0;
    
%     se1 = strel('sphere',2);
%     ctsub = imopen(ctsub,se1);
%     se2 = strel('sphere',5);
%     ctsub = imclose(ctsub,se2);    
    
%     ctsub = imbinarize(ctsub); % leads to MATLAB crash if not assigned on LHS?

    temp_x0 = x; % to compute the lower bounding box horizontal line
    for k=1:z
        for i=1:x
            if(sum(ctsub(i,:,k)) > sum_thresh)
                if i < temp_x0
                    num = num+1;
                    temp_x0 = i;
%                             fprintf('%d, %d\n',temp_x0,k)
                    break
                end
            end
        end
    end
%     x0(q+1) = temp_x0;

    temp_x1 = 0; % to compute the left bounding box vertical line
    arr = zeros(1,z);

    for k=1:z
        for i=x:-1:1
            if(sum(ctsub(i,:,k)) > sum_thresh)
                if i > temp_x1
                    arr(num+1) = temp_x1;
                    num = num+1;
                    temp_x1 = i;
%                             fprintf('%d, %d\n',temp_x1,k)
                    break
                end
            end
        end
    end
%     x1(q+1) = temp_x1;

    temp_y0 = y; % to compute the left bounding box vertical line
    for k=1:z
        for i=1:y
            if(sum(ctsub(:,i,k)) > sum_thresh)
                if i < temp_y0
                    num = num+1;
                    temp_y0 = i;
%                             fprintf('%d, %d\n',temp_y0,k)
                    break
                end
            end
        end
    end
%     y0(q+1) = temp_y0;

    temp_y1 = 0; % to compute the right bounding box vertical line
    for k=1:z
        for i=y:-1:1
            if(sum(ctsub(:,i,k)) > sum_thresh)
                if i > temp_y1
                    num = num+1;
                    temp_y1 = i;
%                             fprintf('%d, %d\n',temp_y1,k)
                    break
                end
            end
        end
    end
%     y1(q+1) = temp_y1;
    ctsub_new = ctsub(temp_x0:temp_x1, temp_y0:temp_y1, :); % cropping

%     fprintf("%d  %d  %d\n", q,q,q);
    fprintf("%d %d %d %d",temp_x0, temp_x1, temp_y0, temp_y1);
    fprintf("********************\n");
    
%     BW = bwconncomp(ctsub_new);
%     BW.NumObjects
%     i = 1;
%     j = 10;
%     k = 1;

%     allprop = regionprops3(BW,'all');
%     surfarea = sort([allprop.SurfaceArea], 'descend');
%     surfmax(q+1) = surfarea(1);
    
    switch (q+1)
        case 1
            save cropped001-100.mat ctsub_new;
        case 2
            save cropped101-200.mat ctsub_new;
        case 3
            save cropped201-300.mat ctsub_new;
        case 4
            save cropped301-400.mat ctsub_new;
        case 5
            save cropped401-500.mat ctsub_new;
        case 6
            save cropped501-600.mat ctsub_new;
        case 7
            save cropped601-700.mat ctsub_new;
        case 8
            save cropped701-800.mat ctsub_new;
        case 9
            save cropped801-900.mat ctsub_new;
        case 10
            save cropped901-1000.mat ctsub_new;
        case 11
            save cropped1001-1100.mat ctsub_new;
        case 12
            save cropped1101-1200.mat ctsub_new;
        case 13
            save cropped1201-1300.mat ctsub_new;
        case 14
            save cropped1301-1400.mat ctsub_new;
        case 15
            save cropped1401-1500.mat ctsub_new;
        case 16
            save cropped1501-1600.mat ctsub_new;
        case 17
            save cropped1601-1700.mat ctsub_new;
        case 18
            save cropped1701-1800.mat ctsub_new;
        case 19
            save cropped1801-1900.mat ctsub_new;
        case 20
            save cropped1901-2000.mat ctsub_new;
        case 21
            save cropped2001-2100.mat ctsub_new;
        case 22
            save cropped2101-2200.mat ctsub_new;
        case 23
            save cropped2201-2300.mat ctsub_new;
        case 24
            save cropped2301-2400.mat ctsub_new;
        case 25
            save cropped2401-2500.mat ctsub_new;
        case 26
            save cropped2501-2600.mat ctsub_new;
        case 27
            save cropped2601-2700.mat ctsub_new;
        case 28
            save cropped2701-2800.mat ctsub_new;
        case 29
            save cropped2801-2900.mat ctsub_new;
        case 30
            save cropped2901-3000.mat ctsub_new;
        case 31
            save cropped3001-3100.mat ctsub_new;
        case 32
            save cropped3101-3200.mat ctsub_new;
        case 33
            save cropped3201-3300.mat ctsub_new;
        case 34
            save cropped3301-3400.mat ctsub_new;
        case 35
            save cropped3401-3500.mat ctsub_new;
        otherwise
            save cropped3501-3600.mat ctsub_new;
    end
    
    clear ctsub;
    clear ctsub_new;
    %toc
end

% x0_sub = max(x0);
% x1_sub = max(x1);
% y0_sub = max(y0);
% y1_sub = max(y1);

toc

%%

tic
ct1 = ct(:,:,1501:1600);

ct1(ct1<12)=0;
ct1(ct1>25)=0;
% ct2 = ct1;

% volshow(ct1)

% figure

% se1 = strel('sphere',2);
% ct1 = imopen(ct1,se1);
% se2 = strel('sphere',5);
% ct1 = imclose(ct1,se2);

% volshow(ct1)


toc;