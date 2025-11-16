
% SAVECODE saves the code into the code folder
% output:
% -------
% This code file is saved into the code folder ("/data/code/").
function [ ] = saveCode()

global sub_num session 

try
    % Get all the matlab files in the directory:
    fileStruct = dir('*.m');

    % Create the save directory:
    directory = fullfile(pwd, 'data',['sub-', num2str(sub_num)],...
        ['ses-',num2str(session)], 'code');

    if ~exist(char(directory),'dir')
        mkdir(directory);
    end
    
    % Copy each file in the subject's folder:
    for i = 1 : length(fileStruct)
        k = 0;
        k = strfind(fileStruct(i).name,'.m');
        if (k ~= 0)
            fileName = fileStruct(i).name;
            source = fullfile(pwd,fileName);
            destination = fullfile(directory,fileName);
            copyfile(source,destination);
        end
    end
    
    % Copy the command line log together with the code:
    logfileName = 'log_recon_time.txt';
    logsource = fullfile(pwd,logfileName);
    logdestination = fullfile(directory,logfileName);
    copyfile(logsource,logdestination);
    
    % Saving the helper functions
    helperFunctionFile = fullfile(pwd,'functions');
    destination = fullfile(directory,'functions');
    copyfile(helperFunctionFile,destination)
    
catch  % Try again if something went wrong:
    fileStruct = dir('*.m');
    
    % Create the save directory:
    directory = fullfile(pwd, 'data',['sub-', num2str(sub_num)],...
        ['ses-',num2str(session)], 'code');
    if ~exist(char(directory),'dir')
        mkdir(directory);
    end
    
    for i = 1 : length(fileStruct)
        k = 0;
        k = strfind(fileStruct(i).name,'.m');
        if (k ~= 0)
            fileName = fileStruct(i).name;
            source = fullfile(pwd,fileName);
            destination = fullfile(directory,fileName);
            copyfile(source,destination);
        end
    end
    
    % Saving the log file:
    logfileName = 'log_wp4_beh.txt';
    logsource = fullfile(pwd,logfileName);
    logdestination = fullfile(directory,logfileName);
    copyfile(logsource,logdestination);
    
    % Saving the helper functions
    helperFunctionFile = fullfile(pwd,'functions');
    destination = fullfile(directory);
    copyfile(helperFunctionFile,destination) % XXX commented out as it crashed
    
end
end