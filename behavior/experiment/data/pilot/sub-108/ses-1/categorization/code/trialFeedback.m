function [blk] = trialFeedback(blk_mat, task_mat)

global spaceBar RestartKey

blk = blk_mat.block(1);

 % generate end of block message
        block_msg = ['End of block ', num2str(blk) ' of ', num2str(task_mat.block(end))];
        next_blk_msg = 'Press space to proceed with experiment';

        if strcmp(blk_mat.task, 'categorization')

            if isnan(blk_mat.RT(end))
                % compute accuracy
                [blk_mat] = compute_performance(blk_mat);
            end

            % get mean accuracy
            mean_acc = mean(blk_mat.trial_accuracy);

            % generate feedback message about accuracy
            accuracy_msg = ['Your accuracy is: ', num2str(round(mean_acc*100)), '%'];

            % give option to repeat practice
            if blk_mat.is_practice(1) && mean_acc < 0.7

                % generate practice feedback message
                next_blk_msg = 'Press R to repeat practice';

            end
        else

            % leave empty if not categorization task
            accuracy_msg = [];
        end 


        % contruct final feedback message
        feedback_msg = [block_msg, newline, newline,...
            accuracy_msg, newline, newline,...
            next_blk_msg];

        % show generated end of block message
        showMessage(feedback_msg);

        % wait for response ("r" = repeats practice, "space bar" = proceeds)
        key_pressed = 0;
        while key_pressed == 0
            [~, ~, Resp] = KbCheck();

            if Resp(RestartKey)
                key_pressed = 1;
                blk = blk - 1;
            elseif Resp(spaceBar)
                key_pressed = 1;
            else
                key_pressed = 0;
            end
        end
end