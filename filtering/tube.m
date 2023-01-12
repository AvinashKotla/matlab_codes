%% Binary 3D model (binarymodel) - single cylinder

height = 512;
width = 512;
nstack = 100
radius = 10;
vol=zeros(height,width,nstack);

vol2 = vol;
slice = vol(:,:,10);
% f = 1.5*(x-x0)^2 + 1.5*(y-y0)^2;
% x0 = randi([0 10],1);
% y0 = randi([0 8],1);

x0 = height/2;
y0 = width/2;

for x=1:height
    for y=1:width
        for z=1:nstack
%             if (x-x0)^2 + (y-y0)^2 < 25*randi([1 4],1)
            if (x-x0)^2 + (y-y0)^2 < radius^2
                vol(x,y,z) = 255;
            end
        end
    end
end

binvol = imbinarize(vol);

volshow(vol)
% volshow(binvol)

% montage({vol(:,:,30),vol(:,:,45),binvol(:,:,30),binvol(:,:,45)});

%% draw circles along a line

x1 = -50:1:-1;
x2 = zeros(1,50);
x3 = 1:1:50;
x = [x1 x2 x3];
y1 = 49:-1:1;
y2 = 0:-1:-50; 
y3 = 0:1:49;
y = [y1 y2 y3];

r = 10;
plot(x,y)
hold on
for i=1:length(y)
    viscircles([x(i),y(i)],r);
%     viscircles([x(i),y(i)],r);
    hold on
end

%% 2d delaunay triangulation
% rng default;
x = rand([20,1]);
y = rand([20,1]);
DT = delaunay(x,y);
triplot(DT,x,y);

%%
x1 = -50:1:-1;
x2 = zeros(1,50);
x3 = 1:1:50;
x = [x1 x2 x3];
y1 = 49:-1:1;
y2 = 0:-1:-50;
y3 = 1:1:50;
y = [y1 y2 y3];


% plot(x,y)

x0 = randi([0 10],1);
y0 = randi([0 8],1);

vol=zeros(length(x),length(x),50);

for i=1:length(x)
    for k=1:50
%         if (x(i)-x0)^2 + (y(i)-y0)^2 > 25*randi([1 4],1)
        if ((x(i).^2) + (y(i).^2)) < 25
            vol(i,i,k) = 255;
        end
    end
end


volshow(vol)

%% PDE solution

alpha = 1;
l = pi;
t = 1:1:10;
n = linspace(0,100);

syms x;
syms t;
p = n;

T = 32/9*pi * (sind(n*x) * exp((-alpha*n^2)*t));

%% Binary 3D model (binarymodel) - branch
tic
height = 512;
width = 512;
nstack = 100;
radius = 10;
vol=zeros(height,width,nstack);
slice = zeros(height,width);

% [x0,y0] = branch(0,0,30,60,200);
x0 = 0:1:100;
y0 = 2*x0+3;
plot(x0,y0)

% x0 = [50 100 150];
% y0 = [50 100 150];

for x=1:height
    for y=1:width
        for z=1:nstack
            for aa=1:length(x0)
                if (x-x0(aa))^2 + (y-y0(aa))^2 < radius^2
                    vol(x,y,z) = 128;
                end
            end
        end
    end
end

% for aa=1:length(x0)
%     for x=1:height
%         for y=1:width
%             for z=1:nstack
%                 if (x-x0(aa))^2 + (y-y0(aa))^2 < radius^2
%                     vol(x,y,z) = 255;
%                 end
%             end
%         end
%     end
% end

volshow(vol)
% imshow(slice)
toc

%%

im = dicomread("chloro_z3/chlorophyta_19_20220514_145135_0150.dcm");
imad = imadjust(im);
imad_copy = imad;
CC = bwconncomp(imad);

numPixels = cellfun(@numel,CC.PixelIdxList);
[smallest,idx] = min(numPixels);
imad(CC.PixelIdxList{idx}) = 0;

montage({imad,imad_copy});

%%
im = dicomread("chloro_z3/chlorophyta_19_20220514_145135_0150.dcm");
imad = imadjust(im);
imad_copy = imad;

% imc = imcontour(imad);
imclevel = imcontour(imad,4);

% montage({imcontour(imad),imcontour(imad,2)})

%%

grid = zeros(100,100);
gridold = grid;
for i=1:100
    for j=1:100
        grid(i,j) = randi([0,2],1);
    end
end

montage({gridold,grid})

CC = bwconncomp(imad);


%% morphological operations

tic

[x0,y0] = branch(0,0,30,60,200);

maxx0 = ceil(max(x0(:))));
minx0 = ceil(min(x0(:))));
maxy0 = ceil(max(y0(:))));
miny0 = ceil(min(y0(:))));
buffer = 100;

height = abs(maxx0)+abs(minx0)+2*buffer;
width = abs(maxy0)+abs(miny0)+2*buffer;

qq = zeros(width,height);
curve = [x0;y0];
[m,n] = size(qq);


% x0 = [200 100 350];
% y0 = [150 50 200];
% radius = [20 30 40];

% for i=1:m
%     for j=1:n
%         for x=1:length(x0)
%             for y=1:length(y0)
%                 if curve(x,y) == [i,j]
%                     qq(i,j) = curve(x,y);
%                 end
%             end
%         end
%     end
% end

qqorig = qq;
qq2 = qq;
qq3 = qq;
qq4 = qq;

% se = strel('line',20,45);
% qq = imdilate(qq,se);
% 
% se2 = strel('disk',20);
% qq2 = imdilate(qq2,se2);
% 
% se3 = strel('octagon',3);
% qq3 = imdilate(qq3,se3);
% 
% montage({qqorig,qq,qq2,qq3});
toc

%%
grid = zeros(512,512,100);
grid(50:55,50:55,20:80) = 50;
grid2 = grid;
se = strel('sphere',100);
grid = imdilate(grid,se);
volshow(grid);


%%
% ni = niftiread('/home/bebm/AKKM/IISER_CT/IISER_Pune_CT/chlorophyta/chloro_STL_ITKSNAP_NII/z2_chloro_NII.nii');
% % volshow(ni);
% 
% reg = regionprops3(ni);
% conn = bwconncomp(ni);
% id = conn.PixelIdxList;

height = 128;
width = 128;
img = zeros(height,width);

img_old = img;
x0 = [1 10];
y0 = [1 20];
r = 3;
for m=1:length(x0)
    for n=1:length(y0)
        for i=1:width
            for j=1:height
                if ((i-x0(m))^2 + (j-y0(n))^2 < r^2)
                    img(i,j) = 128;
                end
            end
        end
    end
end

% montage({img_old,img});
imshow(img)

bw = bwconncomp(img);
ids = bw.PixelIdxList;



    
    
    
    
    
