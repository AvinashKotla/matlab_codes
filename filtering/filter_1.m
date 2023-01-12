I = imread("axial_0427.png");
imgray = imadjust(rgb2gray(I));
imshow(imgray);
% imdistline;
hold on

%% Canny Edge Detection

N = edge(I, 'Canny');
subplot(1, 3, 1),
imshow(N);
title("Canny");

subplot(1,3,2),
imshow(imadjust(I));
title("Adjusted");

subplot(1,3,3),
con = imcontour(I,5);
title("Imcontour");

%% imfindcircles
[centers, radii, metric] = imfindcircles(I,[6 80]);
% centersStrong5 = centers(1:2,:); 
% radiiStrong5 = radii(1:2);
% metricStrong5 = metric(1:2);
% viscircles(centersStrong5, radiiStrong5,'EdgeColor','b');

viscircles(centers, radii);

%% regionprops
stats = regionprops('table',imgray,'Centroid','MajorAxisLength','MinorAxisLength');
% centers = stats.Centroid;
% diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
% radii = diameters/2;
% viscircles(centers,radii*0.5);




