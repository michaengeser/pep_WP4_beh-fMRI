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


% make files names 
if length(category) > 1
    dir = string(fullfile(pwd,'data',['sub-', num2str(sub_num)],...
        ['ses-',num2str(session)], task, category));
    fileName = fullfile(dir, ['sub-', num2str(sub_num), '_ses-', num2str(session),...
        '_run-', num2str(blk_num), '_task-', task, '_', category, '_events']);
else
    dir = string(fullfile(pwd,'data',['sub-', num2str(sub_num)],...
        ['ses-',num2str(session)], task));
    fileName = fullfile(dir, ['sub-', num2str(sub_num), '_ses-', num2str(session),...
        '_run-', num2str(blk_num), '_task-', task, '_events']);
end


% Creating the directories if they don't already exist:
if ~exist(dir, 'dir')
    mkdir(dir);
end

if isnumeric(blk_num)
    if blk_num < 0
        blk_num = 0;
    end
    blk_num = num2str(blk_num);
end


% check if they already exist and add repetition to file name
repetition = 0;
new_fileName = fileName; 

while exist(new_fileName + '.mat', 'file') 
    repetition = repetition+1;
    new_fileName = insertAfter(fileName,'events',['_repetition_',num2str(repetition)]);
end 
fileName = new_fileName;

% save files
save(fileName + '.mat', 'input_table');
writetable(input_table, fileName + '.csv');

end