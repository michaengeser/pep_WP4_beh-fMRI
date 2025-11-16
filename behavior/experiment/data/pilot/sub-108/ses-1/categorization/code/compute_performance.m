function [trial_mat] = compute_performance(trial_mat)
%%

% keys
global KITCHEN_KEY BATHROOM_KEY

%%
disp('WELCOME TO compute_performance');


trial_mat = trial_mat(~isnan(trial_mat.texture), :);

for tr = 1:length(trial_mat.trial)

    trial_mat.trial_accuracy(tr) = nan;
    if strcmp(trial_mat.task{tr}, 'categorization')

        % compute correctness for categorization task task
        if (trial_mat.trial_response(tr) == KITCHEN_KEY && strcmp(trial_mat.category{tr}, 'kitchen')) ||...
                (trial_mat.trial_response(tr) == BATHROOM_KEY && strcmp(trial_mat.category{tr}, 'bathroom'))
            trial_mat.trial_accuracy(tr) = 1;

        elseif (trial_mat.trial_response(tr) == BATHROOM_KEY && strcmp(trial_mat.category{tr}, 'kitchen')) ||...
                (trial_mat.trial_response(tr) == KITCHEN_KEY && strcmp(trial_mat.category{tr}, 'bathroom'))
            trial_mat.trial_accuracy(tr) = 0;
        end
    end


    % get RT 
    trial_mat.RT(tr) = trial_mat.time_of_resp(tr) - trial_mat.stim_time(tr);

end
end
