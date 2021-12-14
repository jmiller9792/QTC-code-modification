function [codeNumString] = getCodeNum(programLine,i)
    splitStr = strsplit(programLine(i+1:end),' ');
    codeNumString = splitStr{1};
end