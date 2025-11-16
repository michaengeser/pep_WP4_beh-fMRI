function [trial_mat] = prepare_log(trial_mat)
% This function adds column to the trial matrix for each event we might
% want to log:

trial_mat.stim_time(:) = nan;  % Time stamp of the visual stimulus onset
trial_mat.time_of_resp(:) = nan;  % Time stamp of response to visual stimulus
trial_mat.trial_response(:) = nan;  % First button pressed
trial_mat.blank_time(:) = nan;  % Time stamp of fixation onset
trial_mat.mask_time(:) = nan;  % Time stamp of fixation onset
trial_mat.fix_time(:) = nan;  % Time stamp of fixation onset
trial_mat.jit_onset(:) = nan;  % Time stamp of the jitter onset
trial_mat.trial_end(:) = nan;  % Time stamp of trial end
trial_mat.trial_accuracy(:) = nan;  % Time stamp of the jitter onset
trial_mat.RT(:) = nan;  % Time stamp of trial end

end

