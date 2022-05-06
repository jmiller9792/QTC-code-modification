
computer = 'JDMLaptop';
computer = 'Purdue';
cd Z:\equip\QTC\Programs

if strcmp(computer,'JDMLaptop')
    addpath('C:\Users\jmill\OneDrive - purdue.edu\Documents\GitHub\QTC-code-modification')
    addpath('C:\Users\jmill\OneDrive - purdue.edu\Documents\GitHub\QTC-code-modification\cnc_edit')
elseif strcmp(computer,'Purdue')
    addpath('U:\My Documents\GitHub\QTC-code-modification')
    addpath('U:\My Documents\GitHub\QTC-code-modification\cnc_edit')
elseif strcmp(computer,'Home')    
    addpath('C:\Users\Justin\Documents\GitHub\QTC-code-modification')
    addpath('C:\Users\Justin\Documents\GitHub\QTC-code-modification\cnc_edit')
end


    