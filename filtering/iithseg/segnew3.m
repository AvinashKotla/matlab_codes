clear;
clc;

%% nifti read
tic
ni = niftiread('/home/biom3d/akkm/medical_imaging/filtering/z2_chloro_NII.nii');
% volshow(ni);

bw = bwconncomp(ni);
prop = regionprops3(bw,'all');
sa = sort([prop.SurfaceArea], 'descend');
samax = sa(1);

toc
%%
tic
v = dicomreadVolume('/home/biom3d/akkm/medical_imaging/filtering/chloro_z3/');
V = squeeze(v);
lower = 2900;
higher = 3980;

V(V<lower) = 0;
V(V>higher) = 0;
V = imbinarize(V);
drawnow; pause(0.05);


BW = bwconncomp(V);
BW.NumObjects
i = 1;
j = 10;
k = 1;

allprop = regionprops3(BW,'all');
surfarea = sort([allprop.SurfaceArea], 'descend');
samaxV = surfarea(1);

% while(and((BW.NumObjects > 5),(samaxV<750000)))
while(BW.NumObjects > 10)
    seopen = strel('sphere',i);
    V = imopen(V,seopen);
    seclose = strel('sphere',j);
    V = imclose(V,seclose);
    
    V = bwareaopen(V,50000);
    
    BW = bwconncomp(V);
    numregions = BW.NumObjects;
    allprop = regionprops3(BW,'all');
    surfarea = sort([allprop.SurfaceArea], 'descend');
    samaxV = surfarea(1);
    i=i+1;
    j=j-1;

    k = k+1;
end
i
j
% V = bwmorph3(V,'clean');
% BW = bwconncomp(V);

% allprop = regionprops3(BW,'all');
% surfarea = sort([allprop.SurfaceArea], 'descend');
% samaxV = surfarea(1);
% numregions = BW.NumObjects


volshow(V)

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
        