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
                coordList = {'X','Y','Z','B','C','I','J','K'};
                for j = 1:length(coordList)
                    if isfield(i.coord,coordList{j}) % Check that coordinate exists in given line
                        if ~isnan(i.coord.(coordList{j}))  % Check that field has valid balue
                            if strcmp(coordList{j},'Z') % Check out of bounds
                                if i.coord.Z>40
                                    i.coord.Z = 40;
                                    i = gCodeCirclesFcn(i,0); % ReCheck circular coordinates
                                    disp([i.lineNum,': Z out of bounds, corrected to 40'])
                                elseif i.coord.Z<-120
                                    i.coord.Z = -120;
                                    i = gCodeCirclesFcn(i,0); % ReCheck circular coordinates
                                    disp([i.lineNum,': Z out of bounds, corrected to -120'])
                                end
                            end
                            lineNew = [lineNew,' ',coordList{j},'=',num2str(i.coord.(coordList{j}),'%.3f')];
                        else
                            disp(['error non-numeric line: ',i.lineNum, ' ', coordList{j}])
                            return
                        end
                    end
                end
%                 if ~isnan(i.coord.X)
%                     lineNew = [lineNew,' X=',num2str(i.coord.(coordList{j}),'%.3f')];
                    
                    
%                 coordList = fieldnames(i.coord);
%                 for j = 1:length(coordList)
%                     if ~isnan(i.coord.(coordList{j})) 
%                         lineNew = [lineNew,' ',coordList{j},'=',num2str(i.coord.(coordList{j}),'%.3f')];
%                     else
%                         disp(['error non-numeric line: ',i.lineNum, ' ', coordList{j}])
%                         return
%                     end
%                     
%                 end

            end
        end
        lineNew = [lineNew,i.tail];
    end
end