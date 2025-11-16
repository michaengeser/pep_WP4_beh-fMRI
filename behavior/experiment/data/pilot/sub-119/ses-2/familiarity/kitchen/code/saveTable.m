% SAVELOG_TABLE saves a matrix as mat and csv
% input:
% ------
% input_table - the table to be saved
% and block numner (can be number or 'all')
%
% output:
% -------
% A mat and a csv file containning the data of log_table are saved to data
% folder


function [] = saveTable(input_table, blk_num)

global sub_num session category task

% Creating the directories if they don't already exist:
dir = string(fullfile(pwd,'data',['sub-', num2str(sub_num)],['ses-',num2str(session)], task, category));
if ~exist(dir, 'dir')
    mkdir(dir);
end

if isnumeric(blk_num)
    if blk_num < 0
        blk_num = 0;
    end
    blk_num = num2str(blk_num);
end

if length(category) > 1
    category_name = ['_', category];
else 
    category_name = [];
end 

% make files names 
fileName = fullfile(dir, ['sub-', num2str(sub_num), '_ses-', num2str(session),...
    '_run-', num2str(blk_num), '_task-', task, category_name, '_events']);

% check if they already exist and add repetition to file name
repetition = 0;
new_fileName = fileName; 

while exist(new_fileName, 'file') 
    repetition = repetition+1;
    new_fileName = insertAfter(fileName,'events',['_repetition_',num2str(repetition)]);
end 
fileName = new_fileName;

% save files
save(fileName + '.mat', 'input_table');
writetable(input_table, fileName + '.csv');

end