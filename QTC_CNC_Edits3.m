% clc
% clear
% close all
function [] = QTC_CNC_Edits3(originalFileName,modifiedFileName,lineNumbering,cConsistency,offsetIndices,offsetValues)
% cd Z:\equip\QTC\Programs
% addpath('U:\My Documents\GitHub\QTC-code-modification')
% addpath('U:\My Documents\GitHub\QTC-code-modification\cnc_edit')
% QTC_CNC_Edits3('nxasbc_57.prg','test.prg',0,0,[],[])


% To be added:
%    Use same function for c consistency as for offset correction?

%Paths for files and functions
% addpath('cnc_edit');
% addpath('Z:\equip\QTC\Programming_Editors\latestQTCprgs');

% %import file to be altered
% originalFileName = 'nxasbc_57.prg';
% modifiedFileName = 'asbc_57.prg'; % doesn't matter which prefix, both will be written

%changes
% lineNumbering = 0;
% cConsistency = 0;
% offsetCorrBool = 0;

prog = fileread(originalFileName);

%% CommentFlip
indexComment = strfind(prog,';');
for i = indexComment
    prog(i) = '(';
end
indexCommentClose = strfind(prog,')');
for i = indexCommentClose
    prog(i) = ' ';
end

%% Split into lines
progLines = strtrim(strsplit(prog,'\n'));
lastCoord = [];
for i = 1:length(progLines)
    temp = parseLine(progLines{i});
    lineStruct(i) = temp;
    if ~isempty(lineStruct(i).coord)
        lineStruct(i).coordLast = lastCoord;
    	lastCoord = lineStruct(i).coord;
    end
end
    
%% Correct Line Numbers (Only uncommented lines. may not even be necessary in QTC)
if lineNumbering
    lineNum = 0;
    for i = 1:length(lineStruct)
        if length(lineStruct(i).lineNum)==4
            lineNumberStr = num2str(lineNum,'%04u');
            lineStruct(i).lineNum = lineNumberStr;
        end
        lineNum = lineNum+10;
    end
end

%% Fix C rotations
% C value can be anything but the machine will move the shortest distance
% to achieve the rotation up to 180deg
if cConsistency
    c_previous = 0; % Start at C = 0deg
    for i = lineStruct
        if isfield(i.coord,'C')
            c_current = i.coord.C;
            c_dif = mod(c_current-c_previous,360); % Remove extra rotations from difference
            if c_dif>180
                c_dif = c_dif-360; % Change range of difference to +-180 deg (may not work??)
            end
            i.coord.C = c_previous+c_dif;
        end
    end
end

%% Offset Correction
% if offsetCorrBool
if ~isempty(offsetIndices)
    for i = 1:length(offsetIndices)
        lineStruct = offsetCorrection(lineStruct,offsetIndices(i),offsetValues(i));
    end
end

%% Create output string\
progMod = [];
for i = lineStruct
    lineNew = writeLine(i);
    if isempty(progMod)
        progMod = char(lineNew);
    else
        progMod = [progMod char(13) char(10) char(lineNew)];
    end
end

%% Apply syntax and filename corrections
header = ['def local path_synch USINT Cnc_X' char(13) char(10)...
    '%001' char(13) char(10)...
    '(SETUP' char(13) char(10)...
    'N0030 G108 ACC=100' char(13) char(10)...
    'N0040 G109 ACC=100' char(13) char(10)];

if strcmp(originalFileName(1:2),'nx') % if original file is from nx, add header.  Will be removed from output nx version
    progMod = [header,progMod];
end
if strcmp(modifiedFileName(1:2),'nx') % if new file name has nx prefix, remove.  This will be added later
    modifiedFileName = modifiedFileName(3:end);
end
%% write new file
fileID = fopen(modifiedFileName,'w');
fprintf(fileID,'%s\n',progMod);
fclose(fileID);

%% Create NX readable version of the file
indexComment = strfind(progMod,'(');
for i = indexComment
    progMod(i) = ';';
end
    
returns = strfind(progMod,10);
progMod = progMod(returns(5)+1:end);
fileID = fopen(['nx',modifiedFileName],'w');
fprintf(fileID,'%s\n',progMod);
fclose(fileID);
end
%% Functions


