% %% Original code - fractal
% function fractal(x1,y1,orientation,levels,i)
% 
% % x2 = x1+cosd(theta)*iterations;
% % y2 = y1+sind(theta)*iterations;
% % i=0;
% x2 = x1+cosd(orientation)*levels;
% y2 = y1+sind(orientation)*levels;
% theta = 30;
% if levels~=0
%     line([x1 x2],[y1 y2],'LineWidth',2);
%     fractal(x2,y2,orientation+theta,levels-1,i);
% %     fprintf("%d %d\n",x2,y2);
%     fprintf("%d\n",3);
%     fractal(x2,y2,orientation-theta,levels-1,i);
%     fprintf("%d",5);
% end
% end


%% fractal_1_theta_input
function [arr_x,arr_y] = fractal(x1,y1,orientation,levels,ltheta,rtheta)

x2 = x1+cosd(orientation)*levels;
y2 = y1+sind(orientation)*levels;

if levels~=0
    line([x1 x2],[y1 y2],'LineWidth',2);
%     fprintf("%d\n",i);
    fractal(x2,y2,orientation+ltheta,levels-1,ltheta,rtheta);
%     fprintf("%d\n",i);
%     fractal(x2,y2,orientation-rtheta,levels-1,ltheta,rtheta,i+1);
end

end

%%
% x=0;
% for i=1:10
%     x = [x,i];
% end
% x
