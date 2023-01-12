%%
clear all;
clc;

arr_x = [];
arr_y = [];

input = {'Enter x value:','Enter y value:', 'Orientation', '#branches', 'leftangle', 'rightangle'};

response = inputdlg(input);

x2 = str2double(response(1));
y2 = str2double(response(2));
orientation = str2double(response(3));
levels = str2double(response(4));
ltheta = str2double(response(5));
rtheta = str2double(response(6));
theta = 30;
% 
% [x,y] = binbranch(0,0,90,5,30,45)
[x,y] = binbranch(x2,y2,orientation,levels);
arr_x = x;
arr_y = y;
% arr_x = [arr_x,x];
% arr_y = [arr_y,y];

levels = levels-1;

% while levels~=0
for i=levels-1:-1:1
    [xl,yl] = binbranch(arr_x,arr_y,orientation+theta,levels);
%     arr_x = [arr_x,xl];
%     arr_y = [arr_y,yl];
    line([arr_x xl],[arr_y yl],'LineWidth',2);
    [xr,yr] = binbranch(xl,yl,orientation-theta,levels);
    arr_x = [arr_x,xr,xl];
    arr_y = [arr_y,yr,yl];
    line([arr_x xr],[arr_y yr],'LineWidth',2);
%     levels = levels-1;
%     
end
% 
% plot(arr_x,arr_y)

