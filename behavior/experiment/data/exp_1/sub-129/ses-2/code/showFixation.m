%SHOWFIXATION - present fixation
% Input: jitter: 1 if you present the fixation for the jitter, 0 if you are
% not in the jitter
% output:
% -------
% fixation_time - the exact time in which it was presented

function [ fixation_time ] = showFixation(Photodiode, feedback_color)


global gray_color w PHOTODIODE white

% if no color is specified its white
if nargin<2
    feedback_color = white;
end

% Draw the fixation
Screen('FillRect', w, gray_color);
drawFixation(feedback_color);
% % Depending on the Photodiode variable, turn the photodiode on or off
if PHOTODIODE
    if strcmp(Photodiode,'PhotodiodeOn')
        drawPhotodiodBlock('on');
        [~,fixation_time] = Screen('Flip', w, [], 1);
    elseif strcmp(Photodiode,'PhotodiodeOff')
        drawPhotodiodBlock('off');
        [~,fixation_time] = Screen('Flip', w, [], 1);
    end
else
    [~,fixation_time] = Screen('Flip', w, [], 1);
end
end
