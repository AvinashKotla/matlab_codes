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

