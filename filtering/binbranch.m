function [x2,y2] = binbranch(x1,y1,orientation,levels)

x2 = x1+cosd(orientation)*levels;
y2 = y1+sind(orientation)*levels;

end