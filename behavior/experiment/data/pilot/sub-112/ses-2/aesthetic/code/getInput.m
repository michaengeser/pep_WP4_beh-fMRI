
% GETINPUT a single iteration of input seeking, returns the RT and address of key pressed
% Input:
% -------
% key - the identity of the key pressed by the user
% Resp_Time - the RT of the user's key press
% PauseTime - If the pause key is pressed, the duration of the pause is
% recorded to make sure we account for it when starting again!
function [ key, Resp_Time ] = getInput()

global NO_KEY ABORT_KEY abortKey KITCHEN_KEY BATHROOM_KEY
global oneKey twoKey threeKey fourKey fiveKey sixKey sevenKey eightKey spaceBar

%%
key = NO_KEY;

[KeyIsDown, Resp_Time, Resp1] = KbCheck();

if KeyIsDown
    if Resp1(oneKey)
        key = 1;
    elseif Resp1(twoKey)
        key = 2;
    elseif Resp1(threeKey)
        key = 3;
    elseif Resp1(fourKey)
        key = 4;
    elseif Resp1(fiveKey)
        key = 5;
    elseif Resp1(sixKey)
        key = 6;
    elseif Resp1(sevenKey)
        key = 7;
    elseif Resp1(eightKey)
        key = 8;
    elseif Resp1(spaceBar)
        key = 9;
    elseif Resp1(KITCHEN_KEY)
        key = KITCHEN_KEY;
    elseif Resp1(BATHROOM_KEY)
        key = BATHROOM_KEY;
    elseif Resp1(abortKey)
        key = ABORT_KEY;
    else
        key = find(Resp1);
    end
end

if key == NO_KEY
    Resp_Time = [];
end
end