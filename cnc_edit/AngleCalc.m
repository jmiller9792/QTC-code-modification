function [ang] = AngleCalc(xtemp,ytemp,xCenter,yCenter)
% angle = atan(y / x) + pi/2;
xtemp = xtemp-xCenter;
ytemp = ytemp-yCenter;

if (xtemp < 0)
    ang = pi - atan(ytemp / -xtemp);
else
    ang = atan(ytemp / xtemp);
end
end