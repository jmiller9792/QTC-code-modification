% clc
% clear
% close all
function [] = QTC_CNC_Edits3(originalFileName,modifiedFileName,lineNumbering,cConsistency,offsetIndices,offsetValues)
% cd Z:\equip\QTC\Programs
% addpath('U:\My Documents\GitHub\QTC-code-modification')
% addpath('U:\My Documents\GitHub\QTC-code-modification\cnc_edit')
% or
% addpath('C:\Users\jmill\OneDrive - purdue.edu\Documents\GitHub\QTC-code-modification')
% addpath('C:\Users\jmill\OneDrive - purdue.edu\Documents\GitHub\QTC-code-modification\cnc_edit')

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
disp('Reading File')
prog = fileread(originalFileName);
disp('File successfully read')
%% CommentFlip
indexComment = strfind(prog,';');
for i = indexComment
    prog(i) = '(';
end
indexCommentClose = strfind(prog,')');
for i = indexCommentClose
    prog(i) = ' ';
end
disp('Comments adjusted')
%% Split into lines
progLines = strtrim(strsplit(prog,'\n'));
% lastCoord = [];
lastCoord.X = 0;
lastCoord.Y = 0;
lastCoord.Z = 0;
circInterpLast = [];
for i = 1:length(progLines)
    temp = parseLine(progLines{i});
    lineStruct(i) = temp;
    
    % if line assigns circular interpolation, store for following lines.
    % If not yet specified, assign xy circular by default
    if ~isempty(lineStruct(i).gNum) % if gNum is specified in a given line
        if strcmp(lineStruct(i).gNum,'17') || strcmp(lineStruct(i).gNum,'18')||...
                strcmp(lineStruct(i).gNum,'19') % check if new circular interpolation
            
            circInterpLast = lineStruct(i).circInterp;
        
        else % if no new ciruclar interpolation
            if isempty(circInterpLast) %and none exists, set default
                circInterpLast = 'xy';
                lineStruct(i).circInterp = circInterpLast;

            else % otherwise take previous
                lineStruct(i).circInterp = circInterpLast;

            end
        end
    %else not included because it interpolation is not necessary for
    %non G codes
    end
    
    if ~isempty(lineStruct(i).coord) % if line has coordinate contents, otherwise do not update last coordinate
        lineStruct(i).coordLast = lastCoord; % Make previous coordinate accessible to current line (for ciruclar interp
%         lineStruct(i).coordLast = lastCoord;
%     	lastCoord = lineStruct(i).coord;
        
        % If coordinate is specified previously but not in current line, add it to current line
        if ~isfield(lineStruct(i).coord,'X')
            lineStruct(i).coord.X = lastCoord.X;
        end
        if ~isfield(lineStruct(i).coord,'Y')
            lineStruct(i).coord.Y = lastCoord.Y;
        end
        if ~isfield(lineStruct(i).coord,'Z')
            lineStruct(i).coord.Z = lastCoord.Z;
        end
        % Save current coordinate to be available for next line
        lastCoord = lineStruct(i).coord;
    end
end
disp('Program split into lines')    
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
disp('Line numbers corrected')
%% Fix C rotations
% C value can be anything but the machine will move the shortest distance
% to achieve the rotation up to 180deg
if cConsistency
    disp('Fixing C Rotations')
    c_previous = 0; % Start at C = 0deg
    for i = 1:length(lineStruct)
        if isfield(lineStruct(i).coord,'C')
            c_current = lineStruct(i).coord.C;
            c_dif = mod(c_current-c_previous,360); % Remove extra rotations from difference
            if c_dif>180
                c_dif = c_dif-360; % Change range of difference to +-180 deg (may not work??)
            end
            lineStruct(i).coord.C = c_previous+c_dif;
            c_previous = lineStruct(i).coord.C;
        end
    end
    
end

%% Offset Correction
% if offsetCorrBool
if ~isempty(offsetIndices)
    disp('Performing offsets')
    for i = 1:length(offsetIndices)
        lineStruct = offsetCorrection(lineStruct,offsetIndices(i),offsetValues(i));
    end
end

%% Create output string\
disp('Generating ouput string')
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
returns = strfind(progMod,10); %carriage returns

header = ['def local path_synch USINT Cnc_X' char(13) char(10)...
    '%001' char(13) char(10)...
    '(SETUP' char(13) char(10)...
    'N0030 G108 ACC=100' char(13) char(10)...
    'N0040 G109 ACC=100' char(13) char(10)];

if strcmp(originalFileName(1:2),'nx') % if original file is from nx, add header.  Will be removed from output nx version
    progMod = [header,progMod];
else % for some reason, there is an error reading the first line which leads to errors writing.  
    progMod = progMod(returns(1)+1:end); % Brute force remove first line
    progMod = ['def local path_synch USINT Cnc_X' char(13) char(10) progMod]; % Brute force add correct first line
end
if strcmp(modifiedFileName(1:2),'nx') % if new file name has nx prefix, remove.  This will be added later
    modifiedFileName = modifiedFileName(3:end);
end
%% write new file
disp('Writing new files')
fileID = fopen(modifiedFileName,'w');
fprintf(fileID,'%s\n',progMod);
fclose(fileID);

%% Create NX readable version of the file
indexComment = strfind(progMod,'(');
for i = indexComment
    progMod(i) = ';';
end
    
returns = strfind(progMod,10); %carriage returns (reset due to changes in header)
progMod = progMod(returns(5)+1:end); % crop from character after fifth carriage return to end
fileID = fopen(['nx',modifiedFileName],'w');
fprintf(fileID,'%s\n',progMod);
fclose(fileID);
end
%% Functions


