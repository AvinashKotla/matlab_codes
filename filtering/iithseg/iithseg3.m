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

zone_size = 100;
zones = ceil(numImages/zone_size);

v = dicomreadVolume('/home/biom3d/akkm/IITH/S4dcm/');
ct = squeeze(v);

toc

%%

tic
% block = 100;
for q=0:30
    diff = (zone_size*q+zone_size)-zone_size*q+1;
    ni2 = zeros(3936,3936,diff);
    [x,y,z] = size(ni2);
    fprintf("%d, %d %d\n",x,y,z);
    fprintf("%d %d\n",zone_size*q+1,zone_size*q+zone_size);
end
toc


%% file cropping

tic
ni2 = zeros(3936,3936,101);
for q=0:12
%     q = 0;
    % ni2 = ct(:,:,zone_size*q+1:zone_size*q+zone_size);
    ni2 = ct(:,:,zone_size*q+1:zone_size*q+zone_size);
    % ni22 = ni2;

    [x,y,z] = size(ni2);
    num=0;
    c=0;
    sum_thresh = 1;

    ni2(ni2 < 13) = 0;
    ni2(ni2 > 25) = 0;
    % imbinarize(ni2);

    temp_x0 = x;
    for k=1:z
    for i=1:x
    %         fprintf("%d\n",i)
        if(sum(ni2(i,:,k)) > sum_thresh)
    %             fprintf("%d, %d, %d\n",i,j,c);
            if i < temp_x0
                num = num+1;
                temp_x0 = i;
    %                     fprintf('%d, %d\n',temp_x0,k)
                break
            end
        end
    end
    end

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
    %                     fprintf('%d, %d\n',temp_x1,k)
                break
            end
        end
    end
    end

    temp_y0 = y;
    for k=1:z
    for i=1:y
    %         fprintf("%d\n",i)
        if(sum(ni2(:,i,k)) > sum_thresh)
    %             fprintf("%d, %d, %d\n",i,j,c);
            if i < temp_y0
                num = num+1;
                temp_y0 = i;
    %                     fprintf('%d, %d\n',temp_y0,k)
                break
            end
        end
    end
    end

    temp_y1 = 0;
    for k=1:z
    for i=y:-1:1
    %         fprintf("%d\n",i)
        if(sum(ni2(:,i,k)) > sum_thresh)
    %             fprintf("%d, %d, %d\n",i,j,c);
            if i > temp_y1
                num = num+1;
                temp_y1 = i;
    %                     fprintf('%d, %d\n',temp_y1,k)
                break
            end
        end
    end
    end

    ni2_new = ni2(temp_x0:temp_x1, temp_y0:temp_y1, :);
    % 
    %     volshow(ni22, 'Colormap',colormap)
    %     figure
    %     hold on
    %     volshow(ni2_new,'Colormap',colormap);

    fprintf("%d  %d  %d\n", q,q,q);
    save cropped101-200.mat ni2_new;


    
%     if (q==1)
%         save cropped1-100.mat ni2_new;
%     else if (q==2)
%           save cropped101-200.mat ni2_new;
%         else
%             save cropped201-300.mat ni2_new;
%         end
%     end
            
    switch q
        case 1
            save cropped1-100.mat ni2_new;
        case 2
            save cropped101-200.mat ni2_new;
        case 3
            save cropped201-300.mat ni2_new;
        case 4
            save cropped301-400.mat ni2_new;
        case 5
            save cropped401-500.mat ni2_new;
        case 6
            save cropped501-600.mat ni2_new;
        case 7
            save cropped601-700.mat ni2_new;
        case 8
            save cropped701-800.mat ni2_new;
        case 9
            save cropped801-900.mat ni2_new;
        case 10
            save cropped901-1000.mat ni2_new;
        case 11
            save cropped1001-1100.mat ni2_new;
        case 12
            save cropped1101-1200.mat ni2_new;
        case 13
            save cropped1201-1300.mat ni2_new;
        otherwise
                save cropped1301-1400.mat ni2_new;
    end
           

%     save ni1.mat ni2_new;

    clear ni2;
    clear ni2_new;

end

toc
