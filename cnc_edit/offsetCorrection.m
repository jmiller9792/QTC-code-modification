function lineStruct = offsetCorrection(lineStruct,offsetKey,offset)
    if ischar(offset) % If offset is specified as character array, convert to number
        offset
        error('ERROR: Offset must be specified as a number')
    end
    for i = 1:length(lineStruct)
        if isfield(lineStruct(i).coord,offsetKey)
            currentValue = lineStruct(i).coord.(offsetKey)
            newValue = currentValue+offset
            lineStruct(i).coord.(offsetKey)=newValue;
        end
        if isfield(lineStruct(i).coordLast,offsetKey)
            currentValueLast = lineStruct(i).coordLast.(offsetKey)
            newValueLast = currentValueLast+offset
            lineStruct(i).coordLast.(offsetKey)=newValueLast;
        end
    end
end