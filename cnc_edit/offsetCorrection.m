function lineStruct = offsetCorrection(lineStruct,offsetKey,offset)
    for i = 1:length(lineStruct)
        if isfield(lineStruct(i).coord,offsetKey)
            currentValue = lineStruct(i).coord.(offsetKey)
            newValue = currentValue+offset
            lineStruct(i).coord.(offsetKey)=newValue;
        end
    end
end