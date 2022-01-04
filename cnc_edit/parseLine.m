function [lineStruct] = parseLine(programLine)
% Converts a character g code line to a structure with pertinent information

% intiialize variables
lineStruct.type = [];
    lineStruct.tail = [];
    lineStruct.lineNum = [];
    lineStruct.gNum = [];
    lineStruct.mCode = [];
    lineStruct.coord = [];
    lineStruct.feed = [];
    lineStruct.coordLast = [];

if isempty(programLine) % Create output
    lineStruct.type = 'blank';
    return
elseif strcmp(programLine(2:4),'def')||strcmp(programLine(1:4),'%001')
    lineStruct.type = 'header';
    lineStruct.tail = programLine;
    return
elseif strcmp(programLine,';def local path_synch USINT Cnc_X')||strcmp(programLine,';%001')
    lineStruct.type = 'header';
    lineStruct.tail = programLine(2:end);
    return
end

% Comment detection and removal
indexComment = strfind(programLine,';');
indexComment = [indexComment,strfind(programLine,'(')];
if ~isempty(indexComment)
    if indexComment(1) == 1
        lineStruct.type = 'comment';
        lineStruct.tail = programLine;
        return
    else
        lineStruct.tail = ' ';
    end
    lineStruct.tail = [lineStruct.tail,programLine(indexComment(1):length(programLine))]; % Commented section returned, only works for line comments, no closure
    programLine = programLine(1:indexComment(1)-1);% Remove from first comment (eliminates issues where comment is not separated by space
end




for i = 1:length(programLine)
    if programLine(i) == 'N'
        lineNum = getCodeNum(programLine,i);
        lineStruct.lineNum = lineNum;
    elseif programLine(i) == 'F'
        feed = getVal2(programLine,i);
        lineStruct.feed = feed;
        lineStruct.type = 'feed';
    elseif programLine(i) == 'G'
        gNum = getCodeNum(programLine,i);
        
        if strcmp(gNum,'100')
            gNum = '01';
        end
        lineStruct.gNum = gNum;
        if strcmp(gNum,'01')
            lineStruct.type = 'linear';
        elseif strcmp(gNum,'02')
            lineStruct.type = 'cwCircle';
        elseif strcmp(gNum,'03')
            lineStruct.type = 'ccwCircle';
        elseif strcmp(gNum,'108')
            lineStruct.type = 'setAccel';
            lineStruct.tail = programLine(i:end);
            return
        elseif strcmp(gNum,'109')
            lineStruct.type = 'setDecel';
            lineStruct.tail = programLine(i:end);
            return
        elseif strcmp(gNum,'17')
            lineStruct.type = 'setInterpPlaneXY';
            lineStruct.tail = [' ',programLine(i:end)];
            return
        elseif strcmp(gNum,'19')
            lineStruct.type = 'setInterpPlaneYZ';
            lineStruct.tail =  [' ',programLine(i:end)];
            return
        else
            disp('G code not recognized');
        end
        for j = i:length(programLine)
            if programLine(j) == 'X' 
                lineStruct.coord.X = getVal2(programLine,j);
            elseif programLine(j) == 'Y'
                lineStruct.coord.Y = getVal2(programLine,j);
            elseif programLine(j) == 'Z'
                lineStruct.coord.Z = getVal2(programLine,j);
            elseif programLine(j) == 'B'
                lineStruct.coord.B = getVal2(programLine,j);
            elseif programLine(j) == 'C'
                lineStruct.coord.C = getVal2(programLine,j);
            elseif programLine(j) == 'I'
                lineStruct.coord.I = getVal2(programLine,j);
                if ~(strcmp(lineStruct.type,'cwCircle')||strcmp(lineStruct.type,'ccwCircle'))
                    disp(['Wrong line type',num2str(lineNum)])
                end                
            elseif programLine(j) == 'J'
                lineStruct.coord.J = getVal2(programLine,j);
                if ~(strcmp(lineStruct.type,'cwCircle')||strcmp(lineStruct.type,'ccwCircle'))
                    disp(['Wrong line type',num2str(lineNum)])
                end       
            elseif programLine(j) == 'K'
                lineStruct.coord.K = getVal2(programLine,j);
            end
        end
    elseif programLine(i) == 'M'
        mCode = getCodeNum(programLine,i);
        lineStruct.type = 'mCode';
        lineStruct.mCode = mCode;
    end
end