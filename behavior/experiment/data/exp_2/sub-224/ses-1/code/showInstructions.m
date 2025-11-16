
%SHOWINSTRUCTIONS - presents the instructions slide
% input:
% ------
% background - the pointer to the instructions image

function [ ] = showInstructions(slide)

    global gray_color w ScreenWidth ScreenHeight PHOTODIODE;
    Screen('FillRect', w, gray_color);

    %show main stimuli
    x = Screen('MakeTexture', w, imread(fullfile(pwd, 'instructions', slide)));
    Screen('DrawTexture',w, x, [], [0 0 ScreenWidth ScreenHeight]);

    if PHOTODIODE
            drawPhotodiodBlock('off')
    end

    Screen('Flip', w);
end
