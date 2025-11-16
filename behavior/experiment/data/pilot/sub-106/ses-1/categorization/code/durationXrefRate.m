% This function makes durations a mutiple of the refresh rate to assure
% optimal timing
%
% Input: trial matrix that controls stimulus presentation
% Output: Updated trial matrix

function [task_mat] = durationXrefRate(task_mat)

global refRate

% loops through all the trials
for tr_time = 1:length(task_mat.trial)

    % duration
    dur_multiplier = round(task_mat.duration(tr_time)/refRate);
    task_mat.duration(tr_time) = refRate*dur_multiplier;

    % blank
    blank_multiplier = round(task_mat.blank(tr_time)/refRate);
    task_mat.blank(tr_time) = refRate*blank_multiplier;
    % mask

    mask_multiplier = round(task_mat.mask_dur(tr_time)/refRate);
    task_mat.mask_dur(tr_time) = refRate*mask_multiplier;
    % jitter

    jit_multiplier = round(task_mat.jitter(tr_time)/refRate);
    task_mat.jitter(tr_time) = refRate*jit_multiplier;

end
end