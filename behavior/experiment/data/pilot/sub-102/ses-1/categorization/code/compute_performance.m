function [blk_mat] = compute_performance(blk_mat)
%%

% keys
global KITCHEN_KEY BATHROOM_KEY

%%
disp('WELCOME TO compute_performance');


blk_mat = blk_mat(~isnan(blk_mat.texture), :);

for tr = 1:length(blk_mat.trial)

    blk_mat.trial_accuracy(tr) = nan;
    if strcmp(blk_mat.task{tr}, 'categorization')

        % compute correctness for categorization task task
        if (blk_mat.trial_response(tr) == KITCHEN_KEY && strcmp(blk_mat.category{tr}, 'kitchen')) ||...
                (blk_mat.trial_response(tr) == BATHROOM_KEY && strcmp(blk_mat.category{tr}, 'bathroom'))
            blk_mat.trial_accuracy(tr) = 1;

        elseif (blk_mat.trial_response(tr) == BATHROOM_KEY && strcmp(blk_mat.category{tr}, 'kitchen')) ||...
                (blk_mat.trial_response(tr) == KITCHEN_KEY && strcmp(blk_mat.category{tr}, 'bathroom'))
            blk_mat.trial_accuracy(tr) = 0;
        end
    end


    % get RT 
    blk_mat.RT(tr) = blk_mat.time_of_resp(tr) - blk_mat.stim_time(tr);

end
end
