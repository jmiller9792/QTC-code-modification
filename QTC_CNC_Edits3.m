% clc
% clear
% close all
function [] = QTC_CNC_Edits3(originalFileName,modifiedFileName,lineNumbering,cConsistency,...
    offsetIndices,offsetValues,feedScale,flipAxis,flipValue,wraps)
% Run set working directories program first

% Flip axis is axis to be flipped about, other is the values of what will
% be flipped.  C_axisZero is detected from first coordinate (set when first
% coordinate is specified).  Flip axis SHOULD NOT BE USED if a different
% calibration check location is used.

% Call function using following command examples
% QTC_CNC_Edits3('nxasbc_57.prg','test.prg',0,0,['X'],[20],1,'',0)
% QTC_CNC_Edits3('nxasbc_57.prg','test.prg',0,0,[],[],1,'',0) 

% cd Z:\equip\QTC\Programs
% addpath('U:\My Documents\GitHub\QTC-code-modification')
% addpath('U:\My Documents\GitHub\QTC-code-modification\cnc_edit')
% or
% addpath('C:\Users\jmill\OneDrive - purdue.edu\Documents\GitHub\QTC-code-modification')
% addpath('C:\Users\jmill\OneDrive - purdue.edu\Documents\GitHub\QTC-code-modification\cnc_edit')

%changes
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
lastCoord.X = 0;
lastCoord.Y = 0;
lastCoord.Z = 0;
circInterpLast = [];
firstCoord = 1;
startWrapIndex = [];
endWrapIndex = [];
finishWrapIndex = [];
for i = 1:length(progLines)
    temp = parseLine(progLines{i});
    lineStruct(i) = temp;
    
    if ~isempty(lineStruct(i).feed)
        lineStruct(i).feed = feedScale*lineStruct(i).feed; % Scale feed
        if lineStruct(i).feed>18000 % implement maximum
            lineStruct(i).lineNum
            lineStruct(i).feed = 18000
        end 
    end

    if strcmp(lineStruct(i).type,'comment')
        if isempty(startWrapIndex)
            if strlength(lineStruct(i).tail)>=6 % check comment length (too short of comments present issues)
                if strcmp(lineStruct(i).tail(2:6),'nWRAP')
                    startWrapIndex = i;
                end
            end
        elseif isempty(endWrapIndex) % if start wrap has already been found, search for end wrap
            if strlength(lineStruct(i).tail)>=9
            if strcmp(lineStruct(i).tail(2:9),'ENDnWRAP')
                endWrapIndex = i;
            end
            end
        else
            if strlength(lineStruct(i).tail)>=8
            if strcmp(lineStruct(i).tail(2:8),'nFINISH')
                finishWrapIndex = i;
            end
            end
        end

    end  
    
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
        if firstCoord
            C_axisZero = lineStruct(i).coord.C
            firstCoord = 0;
        end
        lineStruct(i).coordLast = lastCoord; % Make previous coordinate accessible to current line (for ciruclar interp
        
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
%% mirror program

if strcmp(flipAxis,'X') % swap y values
    error('x axis not implemented yet')
    %C_axisZero should be changed to be parallel to X axis
elseif strcmp(flipAxis,'Y') % swap x values
    for i = 1:length(lineStruct)
        if isfield(lineStruct(i).coord,'X')
            lineStruct(i).coord.X = flipValue-(lineStruct(i).coord.X-flipValue);
            lineStruct(i).coordLast.X = flipValue-(lineStruct(i).coordLast.X-flipValue);
        end
        if isfield(lineStruct(i).coord,'C')
            lineStruct(i).coord.C = C_axisZero-(lineStruct(i).coord.C-C_axisZero);
        end
        if isfield(lineStruct(i).coord,'I')
            lineStruct(i).coord.I = -lineStruct(i).coord.I;
        end
        if strcmp(lineStruct(i).circInterp,'xy')
            if strcmp(lineStruct(i).gNum,'03')
                lineStruct(i).gNum = '02';
            elseif strcmp(lineStruct(i).gNum,'02')
                lineStruct(i).gNum = '03';
            end
        elseif strcmp(lineStruct(i).circInterp,'xz')
            if strcmp(lineStruct(i).gNum,'03')
                lineStruct(i).gNum = '02';
            elseif strcmp(lineStruct(i).gNum,'02')
                lineStruct(i).gNum = '03';
            end
        end
    end
end

%% Wraps
if wraps>1
    % Move finishing code
    for i = finishWrapIndex:length(lineStruct)
        lineStruct((wraps)*(endWrapIndex+1-startWrapIndex)+startWrapIndex+(i-finishWrapIndex))=lineStruct(i);
    end
    
    % Delete extra lines (if previously more than current number of laps)
    startTrimInd = ...
        (wraps)*(endWrapIndex+1-startWrapIndex)+startWrapIndex-1+...
        (length(lineStruct)-(finishWrapIndex-1));
    if startTrimInd<length(lineStruct)
        for i = startTrimInd:length(lineStruct)
            lineStruct(startTrimInd) = [];
        end
    end
    % Copy wrap
    for j = 2:wraps
        % Add numbering to labels
        lineStruct((j-1)*(endWrapIndex+1-startWrapIndex)+startWrapIndex).tail = ...
            strcat('(',lineStruct(startWrapIndex).tail(2:6),num2str(j));
        lineStruct((j-1)*(endWrapIndex+1-startWrapIndex)+endWrapIndex).tail = ...
            strcat('(',lineStruct(endWrapIndex).tail(2:9),num2str(j));
        % Repeat wrap code
        for i = startWrapIndex+1:endWrapIndex
            lineStruct((j-1)*(endWrapIndex+1-startWrapIndex)+i)=lineStruct(i);
        end

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
disp('Line numbers corrected')
%% Fix C rotations
% C value can be anything but the machine will move the shortest distance
% to achieve the rotation up to 180deg
if cConsistency || wraps>1
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
% if isempty(startWrapIndex)
    for i = lineStruct
        lineNew = writeLine(i);
        if isempty(progMod)
            progMod = char(lineNew);
        else
            progMod = [progMod char(13) char(10) char(lineNew)];
        end
    end
% else
%     %Initial movements
%     for i = 1:startWrapIndex-1
%         lineNew = writeLine(lineStruct(i));
%         if isempty(progMod)
%             progMod = char(lineNew);
%         else
%             progMod = [progMod char(13) char(10) char(lineNew)];
%         end
%     end
%     % specified repititions of wraps
%     for j = 1:wraps
%         lineNew = writeLine(lineStruct(startWrapIndex));
%         progMod = [progMod char(13) char(10) char(lineNew) num2str(j)];
%         for i = startWrapIndex+1:endWrapIndex-1
%             lineNew = writeLine(lineStruct(i));
%             progMod = [progMod char(13) char(10) char(lineNew)];
%         end
%         lineNew = writeLine(lineStruct(endWrapIndex));
%         progMod = [progMod char(13) char(10) char(lineNew) num2str(j)];
%     end
%     % Finishing
%     for i = finishWrapIndex:length(lineStruct)
%         lineNew = writeLine(lineStruct(i));
%         progMod = [progMod char(13) char(10) char(lineNew)];
%     end
% end
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
