clc
clear

% To be added:
%    Make more function based.  Allow multiple offset correction in one run
%    Use same function for c consistency as for offset correction
%    Ensure all numbers don't have excessive precision (avoid syntax
%    errors)
%    Remove or comment first lines when comment style is for NX

%Paths for files and functions
addpath('..');
% addpath('gCodeReader'); % getVal.m function used
%import file to be altered
originalFileName = 'asbc_56.prg';
modifiedFileName = 'asbc_56.prg';

%changes
lineNumbering = 0;
cConsistency = 0;
offsetCorrBool = 0;

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
%             lineNum = lineNum+10;
%         elseif isempty(lineStruct(i).lineNum) 
%             continue
%         else
%             disp('Create new case for number field length')
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
if offsetCorrBool
lineStruct = offsetCorrection(lineStruct,'C',-180);
% lineStruct = offsetCorrection(lineStruct,'X',-5)
% lineStruct = offsetCorrection(lineStruct,'Y',71.6--261)
% lineStruct = offsetCorrection(lineStruct,'Z',1.9)
end

%% Create output string\
progMod = [];
for i = lineStruct
    if strcmp(i.type,'header')
        lineNew = i.tail;
    elseif strcmp(i.type,'comment')
        lineNew = i.tail;
    elseif strcmp(i.type,'blank')
        lineNew = [];
    else
        lineNew = strcat('N',i.lineNum);
        if strcmp(i.type,'mCode')
            lineNew = [lineNew,' M',i.mCode];
        elseif strcmp(i.type,'feed')
            lineNew = [lineNew,' F',num2str(i.feed)];
        elseif strcmp(i.type,'setAccel')||strcmp(i.type,'setDecel')
            lineNew = [lineNew,' '];
        elseif strcmp(i.type,'linear')||strcmp(i.type,'cwCircle')||strcmp(i.type,'ccwCircle')
            i = gCodeCirclesFcn(i); % Check circular coordinates
            lineNew = [lineNew,' G',i.gNum];
            if ~isempty(i.coord)
                coordList = fieldnames(i.coord);
                for j = 1:length(coordList)
                    lineNew = [lineNew,' ',coordList{j},'=',num2str(i.coord.(coordList{j}),'%.3f')];
                end
            end
        end
        lineNew = [lineNew,i.tail];
    end
    if isempty(progMod)
        progMod = char(lineNew);
    else
        progMod = [progMod char(13) char(10) char(lineNew)];
    end
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

%% Functions
function lineStruct = offsetCorrection(lineStruct,offsetKey,offset)
    for i = 1:length(lineStruct)
        if isfield(lineStruct(i).coord,offsetKey)
            currentValue = lineStruct(i).coord.(offsetKey)
            newValue = currentValue+offset
            lineStruct(i).coord.(offsetKey)=newValue;
        end
    end
end

