function [tr_mat, file_list] = load_trial_matrix(sub_num)

% open trial matrix
MatFolderName = fullfile(pwd,'trial_matrices');
TableName = ['sub-',num2str(sub_num),'_wp4_beh_trials.csv'];
TablePath = fullfile(MatFolderName, TableName);

% check if this subject has a trial matrix alraedy, if not make one
matrix_folder = [pwd,filesep,'trial_matrices\'];
addpath(matrix_folder)

% make trial matrix
[tr_mat, file_list] = createTrialMatrices(sub_num, 1);

end