clear;
clc;

%%
tic
% v = dicomreadVolume('F:\Work\IISER_Pune_CT\chlorophyta\chloro_highres\chlorophyta_19_20220514_145135');
v = dicomreadVolume(uigetdir);
V = squeeze(v);

input = {'Enter lower threshold:','Enter upper threshold:'};

resp_thresh = inputdlg(input);
lower = str2double(resp_thresh(1));
higher = str2double(resp_thresh(2));

V(V<lower) = 0;
V(V>higher) = 0;
V = imbinarize(V);

input2 = {'Enter the type of morph element:'};
resp_morph = inputdlg(input2);
morph_ele = string(resp_morph(1));

switch morph_ele
    case 'disk'
        input3 = {'Enter the radius'};
        resp_disk = inputdlg(input3);
        rad = str2double(resp_disk(1));
        se = strel(morph_ele,rad);
    case 'diamond'
        input3 = {'Enter the radius'};
        resp_dia = inputdlg(input3);
        rad = str2double(resp_dia(1));
        se = strel(morph_ele,rad);
    case 'square'
        input3 = {'Enter the side'};
        resp_sq = inputdlg(input3);
        sidesq = str2double(resp_sq(1));
        se = strel(morph_ele,sidesq);
    case 'rectangle'
        input3 = {'Enter the height and width'};
        resp_rect = inputdlg(input3);
        h = str2double(resp_rect(1));
        w = str2double(resp_rect(2));
        se = strel(morph_ele,h,w);
    case 'line'
        input3 = {'Enter the length and angle'};
        resp_line = inputdlg(input3);
        len = str2double(resp_line(1));
        deg = str2double(resp_line(2));
        se = strel(morph_ele,len,deg);
    case 'cube'
        input3 = {'Enter the side'};
        resp_cube = inputdlg(input3);
        sidecb = str2double(resp_cube(1));
        se = strel(morph_ele,sidecb);
    case 'sphere'
        input3 = {'Enter the radius'};
        resp_sph = inputdlg(input3);
        radsph = str2double(resp_sph(1));
        se = strel(morph_ele,radsph);
    case 'cuboid'
        input3 = {'Enter the sides'};
        resp_cuboid = inputdlg(input3);
        m = str2double(resp_cuboid(1));
        n = str2double(resp_cuboid(2));
        p = str2double(resp_cuboid(3));
        se = strel(morph_ele,[m n p]);
end

V = imopen(V,se);
volshow(V)

% display the input parameters during volshow

toc

%%
tic
v = dicomreadVolume('F:\Work\IISER_Pune_CT\chlorophyta\chloro_highres\chlorophyta_19_20220514_145135');
V = squeeze(v);
lower = 2900;
higher = 3980;

V(V<lower) = 0;
V(V>higher) = 0;
V = imbinarize(V);
morph_ele = 'disk';
size = 1;

se = strel(morph_ele,size);
se2 = strel('line',35,150);
V = imopen(V,se2);
V = imopen(V,se);
volshow(V)

toc

%%
tic
BW = bwconncomp(V)
allprop = regionprops3(BW,'all');
volsort = sort([allprop.Volume], 'descend');
vlen = length(volsort);
volmax = max(volsort);
volreq = sum(volsort)-sum(volsort(vlen:-1:300));
% V2 = bwareaopen(V,volreq);
se2 = strel('disk',5);
V2 = imerode(V,se2);
volshow(V2)
toc

%%
tic
cube = ones(500,500,500);
x = 1:500;
y = 1:500;
z = 1:500;
r = 20;
% cube(x.^2 + y.^2 + z.^2 < r^2) = 1;
cube(200:300,200:300,200:300) = 255;
toc

volshow(cube)