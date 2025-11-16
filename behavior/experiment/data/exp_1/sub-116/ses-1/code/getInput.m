
% GETINPUT a single iteration of input seeking, returns the RT and address of key pressed
% Input:
% -------
% key - the identity of the key pressed by the user
% Resp_Time - the RT of the user's key press
% PauseTime - If the pause key is pressed, the duration of the pause is
% recorded to make sure we account for it when starting again!
function [ key, Resp_Time ] = getInput()

global NO_KEY valid_resp_keys abortKey kb

%%
key = NO_KEY;

[KeyIsDown, Resp_Time, Resp1] = KbCheck(kb);

if KeyIsDown
    key = find(Resp1);
    key = key(1);


    if key == NO_KEY
        Resp_Time = [];

    elseif key == abortKey
        error('Experiment has been aborted');

        % if pressed key is not a valid response key
    elseif ~min(ismember(key, valid_resp_keys))
        showMessage('Invalid response key!');

        % wait for correct response bar
        [~, ~, wait_resp] = KbCheck(kb);
        while ~wait_resp(valid_resp_keys)
            [~, Resp_Time, wait_resp] = KbCheck(kb);
        end
        key = find(wait_resp);
        key = key(1);
    end
end

end