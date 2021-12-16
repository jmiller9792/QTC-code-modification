function [lineNew] = writeLine(lineStruct)
    i = lineStruct;

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
            i = gCodeCirclesFcn(i,0); % Check circular coordinates
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
end