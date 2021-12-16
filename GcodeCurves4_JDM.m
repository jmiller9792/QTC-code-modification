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

Codes1 = ('N0230 G01 X=-1273.301 Y=268.064 Z=-85.000 B=0.944 C=54.370');
Codes2 = ('N0240 G02 X=-1218.209 Y=310.920 Z=-40.909 B=0.944 C=131.197 I=53.601 J=-12.125 ;r=80');

[line1Struct] = parseLine(Codes1);
[line2Struct] = parseLine(Codes2);
line2Struct.coordLast = line1Struct.coord;
[line2Struct] = gCodeCirclesFcn(line2Struct,1);
disp(writeLine(line2Struct))