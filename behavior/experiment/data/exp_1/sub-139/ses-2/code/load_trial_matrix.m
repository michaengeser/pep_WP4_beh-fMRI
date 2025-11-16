function [tr_mat, file_list] = load_trial_matrix(sub_num)

% open trial matrix
% MatFolderName = fullfile(pwd,'trial_matrices');
% TableName = ['sub-',num2str(sub_num),'_wp4_beh_trials.csv'];
% TablePath = fullfile(MatFolderName, TableName);

% check if this subject has a trial matrix alraedy, if not make one
matrix_folder = [pwd,filesep,'trial_matrices',filesep];
addpath(matrix_folder)

file_name = fullfile(matrix_folder, sprintf("sub-%d_%s_trials.csv", sub_num, 'wp4_beh'));

if exist(file_name,'file')

    % load matrix
    tr_mat = readtable(file_name);

    disp('hello')

    % get file list
    [~, file_list] = createTrialMatrices(sub_num, 1);

else
    % make trial matrix
    [tr_mat, file_list] = createTrialMatrices(sub_num, 1);

     disp('blaaaa')
end

end