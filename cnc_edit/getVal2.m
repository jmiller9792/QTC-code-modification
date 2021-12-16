function [value] = getVal2(line,i)
    % Revised formatting.  (no longer split string)
    % line is a string containing one line of code
    % i is index within the string of the code identifier (x,y,z,...). 
    splitStr = strsplit(line(i+1:end),' '); % Split by spaces after code character
    if isempty(splitStr{1}) % first string is empty if starts with a space
        testString = splitStr{2};
        if strcmp(testString,'=')
            value = str2double(splitStr{3}); %'X = 103.34'
        elseif testString(1)=='='
            value = str2double(testString(2:end)); %'X =103.34'
        else
            disp('error? X 103.34')
        end
    else % if does not start with a space...
        testString = splitStr{1};
        if strcmp(testString,'=')
            value = str2double(splitStr{2});%'X= 103.34'
        elseif strcmp(testString(1),'=')
            value = str2double(testString(2:end)); %'X=103.34'
        else
            if isnumeric(str2double(testString))
                value = str2double(testString); %'X103.34'
            else
                disp('error')
            end
        end
    end
end
  

