%% Circles along a curve
% x = [100,200,300,400];
% y = [100,200,300,400];
r = 2;

y = -50:10:50;
x = sin(y*180/pi);
z = zeros(1,length(y));

plot(x,y,LineWidth=3)
hold on
for i=1:length(y)
    viscircles([x(i),y(i)],r*randi([1,5],1));
%     viscircles([x(i),y(i)],r);
    hold on
end

%%
y = -50:10:50;
x = sin(y*180/pi);

% x = [10,25,30,40];
% y = [100,10,10,50];


slope = zeros(1,length(x)-1);
slope2 = zeros(1,length(x)-1);

for i = 1:length(slope)
    slope(i) = (y(i+1)-y(i))/(x(i+1)-x(i));
    slope2(i) = atan(slope(i));
end
plot(x,y)
hold on
plot(slope)
hold on
plot(slope2,'*')

%%

r=1
teta=-pi:0.01:pi;
x=r*cos(teta);
y=r*sin(teta);
plot3(x,y,zeros(1,numel(x)))

%% Ellipse rotation

a = 2; b = 10;
x = -5:0.1:5;
yminus = -sqrt(b^2 * (1 - x.^2 / a^2));
% plot(x,y)
% hold on
yplus = sqrt(b^2 * (1 - x.^2 / a^2));
plot(x,yminus,x,yplus)

hold on

pplus = zeros(2,length(x));
pplus(1,:) = x;
pplus(2,:) = yplus;

pminus = zeros(2,length(x));
pminus(1,:) = x;
pminus(2,:) = yminus;

theta = 3*pi/6;
rot3d = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
rot2d = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% pnewminus = rot2d * pminus(1:2,:);
% pnewplus = rot2d * pplus(1:2,:);

pnewminus = rot2d * pminus;
pnewplus = rot2d * pplus;

xnewplus = pnewplus(1,:);
ynewplus = pnewplus(2,:);

xnewminus = pnewminus(1,:);
ynewminus = pnewminus(2,:);

plot(xnewplus,ynewplus,xnewminus,pnewminus);

%% Connected regions

slice = zeros(512,512,1);
% x0 = 255;
% y0 = 255;

x0 = [100 200 350 50];
y0 = [120 350 150 400];
r = [30 100 350 650];

for i = 1:length(x0)
    for x=1:512
        for y=1:512
%             if (x-x0(i))^2*randi([-10,10],1) + (y-y0(i))^2  < 400
%             if ((x-x0(i))^2 + (y-y0(i))^2  == 400) && ((x-x0(i))^2 + (y-y0(i))^2  < 400*randi([-100,100],1))
            if ((x-x0(i))^2 + (y-y0(i))^2  < r(i))
                slice(x,y,1) = 255;
            end
        end
    end
end

imshow(slice);

% CC = bwconncomp(slice);
% region = imread("regions2.png");
% r = 5;
% reg = rgb2gray(region);
% imshow(reg);
% 
% hold on
% 
% CC = bwconncomp(reg,4);
% centroid = regionprops(CC,'Centroid');
% all = regionprops(CC,'all');
% 
% 
% for i=1:length(centroid)
%     xc(i) = centroid(i).Centroid(1);
%     yc(i) = centroid(i).Centroid(2);
% end
% 
% % area = all.Area(all.Area > 500);
% 
% for i=1:length(xc)
%     viscircles([xc(i),yc(i)],r);
%     hold on
% end

%% Radon transform main images
% (x0,y0): (0,0) = left top, (512,0) = left bottom, 
% (0,512) = right top, (512,512) = right bottom

slice = zeros(512,512);
% x0 = 255;
% y0 = 255;

x0 = [256 256];
% y0 = [256];
y0 = [128 256];
r = [1600 400]

for i = 1:length(y0)
    for x=1:512
        for y=1:512
            if ((x-x0(i))^2 + (y-y0(i))^2  < r(i))
                slice(x,y) = 255;
            end
        end
    end
end
imwrite(slice,'cir2.png');
imshow(slice);

%%
im1 = imread("cir1.png");
im2 = imread('cir2.png');

p10 = radon(im1,0);
plot(radon(im1,0))
% imshow(iradon(p10,0,'none'),[]) % inverse radon with smearing action
th = 0:45:179;
th2 = 0:2:179;
pangles = radon(im1,th);
% imshow(iradon(pangles,th,'none'),[])

p20 = radon(im2,0);
% imshow(iradon(p20,0,'none'),[]) % smearing action
th = 0:45:179;
pangles2 = radon(im2,th);
% imshow(iradon(pangles2,th,'none'),[])

%% Sinograms

im1 = imread('cir2.png');
im1 = imresize(im1,[100 100]);
iptsetpref('ImshowAxesVisible','on')
theta = 0:180;
[R,xp] = radon(im1,theta);
imshow(R,[],'Xdata',theta,'Ydata',xp,'InitialMagnification','fit');
xlabel('\theta (degrees)')
ylabel('x''')
colormap(gray), colorbar
iptsetpref('ImshowAxesVisible','off');

%%
x = -100:1:100;
y = sin(x);

plot(ff)

%%
rect = zeros(256,256);
rect2 = slice;
rect2(100:150,100:150) = 255;

montage({rect,rect2})

aa = radon(rect2,0);
% plot(radon(rect2,45))
