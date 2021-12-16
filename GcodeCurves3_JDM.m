clc
clear
close all
% This code determines the center of the circle of rotation for a G03 or
% G02 line of G-code. Using I and J values is just one of several ways to
% define an arc.

% A better implementation would allow either pasting of multiple lines
% without needing to assign specific lines of code.  AND/OR code should read
% and entire file and check for issues (any changes to radius would still 
% need to be specified 

% Rotation Direction: G02 = CW,  G03 = CCW
codeDir = 'G02';
r = 55;  %[mm]

% Codes1 = input('line1: ','s');
% Codes2 = input('line2: ','s');
% codeDir = input('codeDir: ','s');
% r = input('codeDir: ','s');

Codes1 = ('N0230 G01 X=-1273.301 Y=268.064 Z=-85.000 B=0.944 C=54.370');
Codes2 = ('N0240 G02 X=-1218.209 Y=310.920 Z=-40.909 B=0.944 C=131.197 I=53.601 J=-12.125 ;r=55');
%  
% (actuate top rear)
% N0730 M56
%N0230 G01 X=-1273.301 Y=268.064 Z=-85.000 B=0.944 C=54.370
%N0240 G01 X=-1218.209 Y=313.920 Z=-50.000 B=0.944 C=131.197 ;r=55
% 
% 

[lineStruct] = parseLine(Codes1)
[lineStruct] = parseLine(Codes2)

Codes1 = strsplit(Codes1,' ');
Codes2 = strsplit(Codes2,' ');
for i = 1: length(Codes1)
    CodeTemp = Codes1{i};
    CodeIndicator = CodeTemp(1);
    if CodeIndicator == 'X'
        x1 = getVal2(Codes1,i);
    elseif CodeIndicator == 'Y'
        y1 = getVal2(Codes1,i);
    end
end
for i = 1: length(Codes2)
    CodeTemp = Codes2{i};
    CodeIndicator = CodeTemp(1);
    if CodeIndicator == 'X'
        x2 = getVal2(Codes2,i);
    elseif CodeIndicator == 'Y'
        y2 = getVal2(Codes2,i);
    end
end
    

% Radius (make sure it is big enough to span the two points)
d = sqrt((x1-x2)^2+(y1-y2)^2);
if (2*r^2) < d
    fprintf('\n   You done messed up. Make the radius bigger\n')
end

x_m = (x2+x1)/2;
y_m = (y2+y1)/2;
m_mid = -(x2-x1)/(y2-y1); %perpendicular to chord between start and end points

b_mid = y_m-m_mid*x_m;
% ymcheck = m_mid*x_m+b_mid

syms xR yR xc yc
test(1) = r^2 == (xc-x1)^2 + (yc-y1)^2;
test(2) = yc == m_mid*xc + b_mid;

T = solve(test,[xc;yc]);

xC1 = double(T.xc(1));
yC1 = double(T.yc(1));
xC2 = double(T.xc(2));
yC2 = double(T.yc(2));


[theta_center1,~] = cart2pol(xC1-x1,yC1-y1);
[theta_center2,~] = cart2pol(xC2-x1,yC2-y1);
[theta2,~] = cart2pol(x2-x1,y2-y1);

thetadif1=theta_center1-theta2;
thetadif2=theta_center2-theta2;

%Correct for out of bounds values
if thetadif1>=pi/2
    thetadif1 = thetadif1-2*pi;
elseif thetadif1<=-pi/2
    thetadif1 = thetadif1+2*pi;
end

if thetadif2>=pi/2
    thetadif2 = thetadif2-2*pi;
elseif thetadif2<=-pi/2
    thetadif2 = thetadif2+2*pi;
end

% Determine which solution corresponds to which direction
if thetadif1>=0
    xC_g03 = xC1;
    yC_g03 = yC1;
    xC_g02 = xC2;
    yC_g02 = yC2;
else
    xC_g02 = xC1;
    yC_g02 = yC1;
    xC_g03 = xC2;
    yC_g03 = yC2;    
end

% Assign 
if strcmp(codeDir,'G02')
    xC = xC_g02;
    yC = yC_g02;
else
    xC = xC_g03;
    yC = yC_g03;
end

I = xC-x1;
J = yC-y1;
% disp(strcat('X',x1,' Y',y1,' I',I,' J',J))
disp(codeDir+ " X" + num2str(x2) + " Y" + num2str(y2) + " I" + num2str(I) + " J" + num2str(J)+ " (R="+ num2str(r)+ ")");

figure()
h1 = axes;
set(h1, 'Ydir', 'reverse')
set(h1, 'Xdir', 'reverse')
hold on
xlabel('X axis, head movement relative to table');
ylabel('Y axis, head movement relative to table');
plot(x1,y1,'b+')
plot(x2,y2,'r+')
plot(xC,yC,'g^')

radius = r;
angStart = AngleCalc(x1,y1,xC,yC);
angEnd = AngleCalc(x2,y2,xC,yC);
if strcmp(codeDir,'G02')
    if angEnd>angStart
        angEnd = angEnd-2*pi;
    end
else
    if angEnd<angStart
        angEnd = angEnd+2*pi;
    end
end

circr = @(radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)];         % Circle Function For Angles In Radians
circd = @(radius,deg_ang)  [radius*cosd(deg_ang);  radius*sind(deg_ang)];       % Circle Function For Angles In Degrees
N = 25;                                                         % Number Of Points In Complete Circle
rad_angl = linspace(angStart,angEnd, N);                             % Angle Defining Arc Segment (radians)
xy_r = circr(radius,rad_angl);                                    % Matrix (2xN) Of (x,y) Coordinates
% figure(1)

hold on
plot(xy_r(1,:)+xC, xy_r(2,:)+yC)                                % Draw An Arc Of Blue Stars
% axis([-1.25*radius  1.25*radius    0  1.25*radius])             % Set Axis Limits
axis equal                                                      % No Distortion With ‘axis equal’

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
