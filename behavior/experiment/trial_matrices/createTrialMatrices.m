% function to create trial matrices
%
% Input the number of trials matrices you want to creat (n) and the first
% subject number (first_sub_num).
%
% The task and categories can be adjusted in first few lines

function [final_mat, file_list] = createTrialMatrices(first_sub_num, n, exp_num)

if nargin < 3
    exp_num = input('Experiment number: [Default: 1]');
    if isempty(exp_num); exp_num = 1; end % default
end

global FILE_POSTFIX expDir

if isempty(FILE_POSTFIX)
    FILE_POSTFIX = '*tif';
end


%% Set constant parameters:
n_pics_per_cate = 50;
n_pics_practice = 10;

% Evaluate experiment categories, tasks, and presentation time
if exp_num == 1
    cate1 = 'bathroom';
    cate2 = 'kitchen';
    tasks = ["categorization", "typicality", "familiarity", "aesthetic", "usability", "complexity"];
    num_frames_stim = 4;
    num_frames_mask = 4;
elseif exp_num == 2
    cate1 = 'bedroom';
    cate2 = 'livingroom';
    tasks = "categorization";
    num_frames_stim = 6;
    num_frames_mask = 4;
end 

pre_file_list = struct(cate1, nan(1,n_pics_per_cate), ...
    cate2, nan(1,n_pics_per_cate), "practice_exp" + num2str(exp_num), nan(1,n_pics_practice));
file_list = struct("SHINEd", pre_file_list, "raw", pre_file_list);
cateories = [cate1, cate2, ["practice_exp" + num2str(exp_num)]];
img_types = ["SHINEd", "raw"];

% Get the list of stimuli files
if isempty(expDir)
    currentDirectory = pwd;
    [expDir, ~, ~] = fileparts(currentDirectory);
end

% initialize loop varaibles
textures = (1:1000) + 10;
count = 0;

for type = img_types
    for cate = cateories
        stimuli_path = fullfile(expDir, 'stimuli', type, cate);
        file_list.(type).(cate) = dir(fullfile(stimuli_path, FILE_POSTFIX));
        for t = 1:length(file_list.(type).(cate))
            count = count + 1;
            file_list.(type).(cate)(t).texture = textures(count);
        end
    end
end

% Define trial matrix table
columnNames = ["sub_num", "trial", "block", "task", "is_practice", "category",...
    "image", "texture", "duration", "mask_dur", "jitter"];
dataTypes = {'double', 'double', 'double', 'string', 'double', 'string',...
    'string', 'double', 'double', 'double', 'double'};

% Create an empty table with the specified column names and data types
final_mat = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);

% Create the jitter distribution:
jitter_mean = 1;
jitter_min = 0.7;
jitter_max = 2;
exp_dist = makedist("Exponential", "mu", jitter_mean);
jitter_distribution = truncate(exp_dist, jitter_min, jitter_max);

%% Creating the trials table:
for sub_num = first_sub_num:first_sub_num + (n-1)

    for task = tasks

        % Set parameters according to task (durations in sec)
        if strcmp(task, 'categorization')

            rng(sub_num) % different random number for each subject
            repeats = 12; % how often each image is shown
            pics_per_blk = 50; % how many images are one block
            duration = (1/60) * num_frames_stim; % for how long is the image shown
            img_type = 'SHINEd'; % images properties are normalized with SHINE toolbox
            mask_dur = (1/60) * num_frames_mask; % mask duration

        else % for rating task

            rng(find(strcmp(task, tasks))) % same random number for each subject but different for each task
            repeats = 1; % how often each image is shown
            pics_per_blk = 25; % how many images are one block
            duration = 1.5; % for how long is the image shown
            img_type = 'raw'; % gray scale images 
            mask_dur = 0; % mask duration

        end

        % Create an empty table with the specified column names and data types
        task_mat = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);

        % get all image names and textures
        count = 0;
        for cate = cateories
            % supress warnings (it will complain about the assign values in the
            % table)
            warning('off','all')
            for t = 1:length(file_list.(img_type).(cate))
                count = count + 1;
                task_mat.sub_num(count) = sub_num;
                task_mat.task(count) = task;
                task_mat.image{count} =  file_list.(img_type).(cate)(t).name;
                task_mat.texture(count) =  file_list.(img_type).(cate)(t).texture;
                task_mat.duration(count) =  duration;
                task_mat.mask_dur(count) =  mask_dur;

                if contains(file_list.(img_type).(cate)(t).name, 'practice')
                    task_mat.is_practice(count) = 1;

                    if contains(task_mat.image{count}, cate2)
                        task_mat.category(count) = cate2;
                    elseif contains(task_mat.image{count}, cate1)
                        task_mat.category(count) = cate1;
                    else
                        task_mat.category(count) = 'unkown';
                    end

                else
                    task_mat.category(count) = cate;
                    task_mat.is_practice(count) = 0;
                end
            end
            % turn back on warnings
            warning('on','all')
        end

        % split task and practice matrix
        practice_mat = task_mat(task_mat.is_practice == 1, :);
        task_mat = task_mat(task_mat.is_practice ~= 1, :);

        for rep = 1:repeats

            % shuffle row order
            task_mat = task_mat(randperm(height(task_mat)),:);

            if rep == 1
                tr_mat = task_mat;
            else
                tr_mat = vertcat(tr_mat, task_mat);
            end
        end

        % shuffle practice
        practice_mat = practice_mat(randperm(height(practice_mat)),:);

        % sort by category if rating task then add block and trial number
        if ~strcmp(task, 'categorization')
            tr_mat = sortrows(tr_mat,'category');
            practice_mat = sortrows(practice_mat,'category');

            blks = repelem(1:(height(tr_mat)/pics_per_blk)/2, pics_per_blk);
            tr_mat.block = repmat(blks, 1, 2)';

        else 
            tr_mat.block = repelem(1:height(tr_mat)/pics_per_blk, pics_per_blk)';
        end
   
        tr_mat.trial = repmat(1: pics_per_blk, 1, height(tr_mat)/pics_per_blk)';

        % add block and trial numbers for practice
        practice_mat.block = zeros(height(practice_mat), 1);
        practice_mat.trial = (1:height(practice_mat))';

        % add new task matrix to trial matrix
        final_mat = vertcat(final_mat, practice_mat);
        final_mat = vertcat(final_mat, tr_mat);

    end

    % Add the jitter:
    final_mat.jitter = random(jitter_distribution, height(final_mat), 1);

    % Save to file:
    % Create file name:
    target_dir = fullfile(expDir, 'trial_matrices');


    file_name = fullfile(target_dir, sprintf("sub-%d_%s_trials.csv", sub_num, 'wp4_beh'));
    writetable(final_mat, file_name);

end

