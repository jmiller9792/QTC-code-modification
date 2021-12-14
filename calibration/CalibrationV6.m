clc
clear
close all
addpath('circFit3d')
addpath('CalibrationPoints')
% prog = fileread('Z:\sbc\Manufacturing\oldPrograms\asbc_mod47_8.prg');

%% Manual data entry (from QTC program)
cValNominal = 0;%107.568; %C or R from the QTC code when T2 tilt is being tested
bValNominal = 0;%B or T2 from the QTC code when T1 tilt is being tested. 
% Positive or negative should be checked...
ind = 2; % Point of T1 motion from which zero is referenced. Use point not near zero to avoid direction issues
t1ValNominal = -19.07; %from the QTC manually when T1 tilt is being tested AT INDEX ind. not available from QTC program

tableRefIndex = 19;
headRefIndex = 20;

%% Data Read
% Generate a file from the leica scanner moving the axes in x,y,z,R,T2,and
% T1. At least 2 line points and 3 rotation points are required
% points = readmatrix('CalibrationPoints.xlsx');
points = readmatrix('calFull4.csv');
pointNum = points(:,1);
h = points(:,2);
v = points(:,3);
d = points(:,4);
%Convert h, v, d to x,y,z coordinates
x = d.*sind(v).*cosd(h); 
y = d.*sind(v).*sind(h);
z = d.*cosd(v);

%% Get XYZ directions
origin = [x(2),y(2),z(2)]; % origin in csys of tracker
XmotionInd = [4,6:9];
YmotionInd = [4,10:14];
ZmotionInd = [4,15:18];
[dirX] = getEigDir(XmotionInd,x,y,z,'L');
[dirY] = getEigDir(YmotionInd,x,y,z,'L');
[dirZ] = getEigDir(ZmotionInd,x,y,z,'L');

T2motionInd = [20:24];
RmotionInd = [24:28];
T1motionInd = [29:33];
[dirR] = getEigDir(RmotionInd,x,y,z,'C');
[dirT2] = getEigDir(T2motionInd,x,y,z,'C');
[dirT1] = getEigDir(T1motionInd,x,y,z,'C');

% XmotionInd = 4:5;
% YmotionInd = 2:3;
% [dirX] = getEigDir(XmotionInd,x,y,z,'L');
% [dirY] = getEigDir(YmotionInd,x,y,z,'L');
% dirZ = cross(dirX.directionVector,dirX.directionVector);

XtableInd = 1:2;
YtableInd = 2:3;
[tabX] = getEigDir(XtableInd,x,y,z,'L');
[tabY] = getEigDir(YtableInd,x,y,z,'L');

t = (-100:10:100);

figure(1)
hold on
plot3(x,y,z); % Plot traced path

plotFitLine(t,dirX,'r')
plotFitLine(t,dirY,'g')
plotFitLine(t,dirZ,'y')

plotCircle3D(t,dirT2,'r')
plotCircle3D(t,dirT1,'g')
plotCircle3D(t,dirR,'b')

axis equal

%% Convert points to new CSYS (_m for machine motion)
%https://math.stackexchange.com/questions/2306319/transforming-point-between-euclidean-coordinate-systems
for i = 1:length(x)
    x_m(i) = ([x(i),y(i),z(i)]-origin)*dirX.directionVector;
    y_m(i) = ([x(i),y(i),z(i)]-origin)*dirY.directionVector;
    z_m(i) = ([x(i),y(i),z(i)]-origin)*dirZ.directionVector;
end

%% T1 zero calc
% Get directiona nd motion centers of T1 and T2
syms xsym ysym zsym
A = dirT2.directionVector(1);
B = dirT2.directionVector(2);
C = dirT2.directionVector(3);
a = dirT2.motionCtr(1);
b = dirT2.motionCtr(2);
c = dirT2.motionCtr(3);

xd = dirT1.directionVector(1);
yd = dirT1.directionVector(2);
zd = dirT1.directionVector(3);
xc = dirT1.motionCtr(1);
yc = dirT1.motionCtr(2);
zc = dirT1.motionCtr(3);
r = dirT1.radius;

% Find points where the T1 arc circle is "lowest" along axis of rotation
% This is accomplished by solving a system of equations consisting of the
% (t1 sphere, t1 plane -> t1 arc) and the plane intersecting the 
eqnT2Plane = A*(xsym-a)+B*(ysym-b)+C*(zsym-c) == 0;
eqnT1Plane = xd*(xc-xsym)+yd*(yc-ysym)+zd*(zc-zsym)== 0;
eqnT1Sphere = (xsym-xc)^2+(ysym-yc)^2+(zsym-zc)^2 == r^2;
eqns = [eqnT2Plane,eqnT1Plane,eqnT1Sphere];
s = solve(eqns,[xsym,ysym,zsym]);
hold on
plot3(double(s.xsym(1,1)),double(s.ysym(1,1)),double(s.zsym(1,1)),'x')
plot3(double(s.xsym(2,1)),double(s.ysym(2,1)),double(s.zsym(2,1)),'v')

if s.zsym(1,1)<s.zsym(2,1)
    zeroPoint = [s.xsym(1,1),s.ysym(1,1),s.zsym(1,1)];
else
    zeroPoint = [s.xsym(2,1),s.ysym(2,1),s.zsym(2,1)];
end
zpVect = double(zeroPoint-dirT1.motionCtr);
actPtVect = double([x(T1motionInd(ind)),y(T1motionInd(ind)),z(T1motionInd(ind))]-dirT1.motionCtr);

%Compute angle between reference and zero point
u = zpVect;
v = actPtVect;
CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
AngleMeasT1 = real(acosd(CosTheta));
t1ZeroOffset = t1ValNominal+AngleMeasT1


%% Offset Measurements
% C offset measurement
AngleMeasC = acosd(dot(dirT2.directionVector,dirX.directionVector));% Angle between T2 axis and X direction
%Initial value used in NX model for C axis
cZeroOffset = 90+cValNominal-AngleMeasC % Add 90 due to current NX model

%% B offset measurement
AngleMeasB = 90-acosd(dot(dirT1.directionVector,dirR.directionVector));% Angle between T1 axis and line perpendicular to R axis
%Initial value used in NX model for B axis
bZeroOffset = bValNominal-AngleMeasB


%% Table
tableRefOffset = 10; % 10mm offset inboard for marker location
tableRefIndexSW = 2;

jigRefIndex = 4;
JigRefX = x_m(jigRefIndex)-x_m(tableRefIndexSW)+tableRefOffset;
JigRefY = y_m(jigRefIndex)-y_m(tableRefIndexSW)+tableRefOffset;
JigRefZ = z_m(jigRefIndex)-z_m(tableRefIndexSW);
JigRef1Coords = [JigRefX,JigRefY,JigRefZ]

jigRefIndex = 5;
JigRefX = x_m(jigRefIndex)-x_m(tableRefIndexSW)+tableRefOffset;
JigRefY = y_m(jigRefIndex)-y_m(tableRefIndexSW)+tableRefOffset;
JigRefZ = z_m(jigRefIndex)-z_m(tableRefIndexSW);
JigRef2Coords = [JigRefX,JigRefY,JigRefZ]


%% XYZ offset
xoffset = x_m(headRefIndex)-x_m(tableRefIndex);
yoffset = y_m(headRefIndex)-y_m(tableRefIndex);
zoffset = z_m(headRefIndex)-z_m(tableRefIndex);
offsetCoords = [xoffset,yoffset,zoffset]
plot3(x(headRefIndex),y(headRefIndex),z(headRefIndex),'v')
plot3(x(tableRefIndex),y(tableRefIndex),z(tableRefIndex),'v')

figure(2)
plot3(x_m,y_m,z_m);
axis equal


%% R axis angles
zenith =  acosd(dot(dirT2.directionVector,dirX.directionVector));
% azimuth = 
% RangleX = acosd(dot(dirX.directionVector,dirR.directionVector))-90 %-90 adjusts to x direction being inverted due to measurement point order
% RangleY = 90-acosd(dot(dirY.directionVector,dirR.directionVector)) %90- adjusts to reference direction in CAD model
RangleX = acosd(dot(dirX.directionVector,dirR.directionVector)); %-90 adjusts to x direction being inverted due to measurement point order
RangleY = acosd(dot(dirY.directionVector,dirR.directionVector)); %90- adjusts to reference direction in CAD model
Rangles = [RangleX,RangleY]

%% Functions

function [dir] = getEigDir(motionInd,x,y,z,type)
% Using the SVD technique to determine the eigenvector of the line of
% points (indexes specified by motionInd, all points in xyz). The highest
% sigma value indicates the best fit of a line('L').  lowest value is normal to
% a circular path ('C')
% https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
    motionPts = [x(motionInd),y(motionInd),z(motionInd)];
    motionCtr = mean(motionPts,1);
    A = motionPts-motionCtr;
    [U,S,V] = svd(A);
    if strcmp(type,'L')
        direction = V(:,1); % line direction
        radius = 0;
    elseif strcmp(type,'C')
        [motionCtr, direction, radius] = CircFit3D(motionPts);
    end
    dir.directionVector = direction;
    dir.motionCtr = motionCtr;
    dir.radius = radius;
end
function plotFitLine(t,dir,color)
    direction = dir.directionVector;
    center = dir.motionCtr;
    hold on
    coord = direction*t+center'*ones(1,length(t));
    plot3(coord(1,:),coord(2,:),coord(3,:),color,'linewidth',1.5)
end
% 
function plotCircle3D(t,dir,color)
    normal = dir.directionVector';
    center = dir.motionCtr;
    radius = dir.radius;
    % https://www.mathworks.com/matlabcentral/fileexchange/26588-plot-circle-in-3d
    theta=0:0.01:2*pi;
    v=null(normal);
    points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    plot3(points(1,:),points(2,:),points(3,:),color); %Plot circle
    plotFitLine(t,dir,color)
    plot3(center(1),center(2),center(3),color);
end