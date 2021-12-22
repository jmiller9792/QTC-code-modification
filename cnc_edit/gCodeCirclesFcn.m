function [line] = gCodeCirclesFcn(line,plotBool)

errorAccept = 0.01; 
% % 
% This code determines the center of the circle of rotation for a G03 or
% G02 line of G-code. Using I and J values is just one of several ways to
% define an arc.

% Rotation Direction: G02 = CW,  G03 = CCW

if strcmp(line.type,'cwCircle')
    codeDir = 'G02';
elseif strcmp(line.type,'ccwCircle')
    codeDir = 'G03';
else
    return
end

x1 = line.coordLast.X;
y1 = line.coordLast.Y;
x2 = line.coord.X;
y2 = line.coord.Y;

% Get radius from line (this is used for quick changes to radius)
comment = line.tail(2:end);
rIndex = strfind(comment,'r');
rIndex = [rIndex,strfind(comment,'R')];
if ~isempty(rIndex)
    r = getVal2(comment,rIndex);
    
else
    disp([line.lineNum,'No radius found'])
    return
end

if ~isempty(line.coord.I) % If already specified, check that radius is consistent
    % Find center coordinates
    xcCurrent = x1+line.coord.I;
    ycCurrent = y1+line.coord.J;
    rad2start = sqrt((x1-xcCurrent)^2+(y1-ycCurrent)^2);
    rad2end = sqrt((xcCurrent-x2)^2+(ycCurrent-y2)^2);
    % if center point is not equal distant from points
      
    if abs(rad2start-rad2end) > errorAccept% || abs(rad2end-r) > errorAccept
        % If unacceptable error, display actual radii
        disp('Radii do not match')
        disp(['Radius to start:' num2str(rad2start)])
        disp(['Radius to end:' num2str(rad2end)])
        disp(['Difference:' num2str(abs(rad2start-rad2end))])
        disp([line.lineNum,' corrected radius collected from comment: ',num2str(r)])
        plotBool = 1;
    else % if consistent no action necessary
        return
    end
    
end

% Radius (make sure it is big enough to span the two points)
d = sqrt((x1-x2)^2+(y1-y2)^2);
if (2*r) <= d
   disp('You done messed up. Make the radius bigger')
   return
end

x_m = (x2+x1)/2;
y_m = (y2+y1)/2;
m_mid = -(x2-x1)/(y2-y1); %perpendicular to chord between start and end points

b_mid = y_m-m_mid*x_m;
% ymcheck = m_mid*x_m+b_mid

syms xc yc
test = r-sqrt((xc-x1)^2 + ((m_mid*xc + b_mid)-y1)^2);
T = solve(test,xc);
xC1 = double(T(1));
xC2 = double(T(2));
yC1 = m_mid*xC1 + b_mid;
yC2 = m_mid*xC2 + b_mid;

% test(1) = r == sqrt((xc-x1)^2 + (yc-y1)^2);
% test(2) = yc == m_mid*xc + b_mid;
% circr = @(xc)  r-sqrt((xc-x1)^2 + ((m_mid*xc + b_mid)-y1)^2)

% T = solve(test,[xc;yc])

% xC1 = double(T.xc(1));
% yC1 = double(T.yc(1));
% xC2 = double(T.xc(2));
% yC2 = double(T.yc(2));


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
% disp(codeDir+ " X" + num2str(x2) + " Y" + num2str(y2) + " I" + num2str(I) + " J" + num2str(J)+ " (R="+ num2str(r)+ ")");
line.coord.I = I;
line.coord.J = J;
disp(['Solved error = ',num2str(abs(sqrt((x1-xC)^2+(y1-yC)^2)-sqrt((xC-x2)^2+(yC-y2)^2)))])
disp([num2str(yC),'=',num2str(m_mid*xC + b_mid)])
disp([num2str(r),'=',num2str(sqrt((xC-x1)^2 + (yC-y1)^2))])

if plotBool
    

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

% plotting
circr = @(radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)];         % Circle Function For Angles In Radians
circd = @(radius,deg_ang)  [radius*cosd(deg_ang);  radius*sind(deg_ang)];       % Circle Function For Angles In Degrees
N = 25;                                                         % Number Of Points In Complete Circle
rad_angl = linspace(angStart,angEnd, N);                             % Angle Defining Arc Segment (radians)
xy_r = circr(radius,rad_angl);                                    % Matrix (2xN) Of (x,y) Coordinates
% figure(1)

hold on
plot(xy_r(1,:)+xC, xy_r(2,:)+yC)                                % Draw An Arc Of Blue Stars
% axis([-1.25*radius  1.25*radius    0  1.25*radius])             % Set Axis Limits
axis equal  % No Distortion With ‘axis equal’
plot(xcCurrent,ycCurrent,'>')
% plot(xC,yC,'x')
legend('Start','End','Center','Path','PreviousCenter')
end
end