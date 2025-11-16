
% GETINPUT a single iteration of input seeking, returns the RT and address of key pressed
% Input:
% -------
% key - the identity of the key pressed by the user
% Resp_Time - the RT of the user's key press
% PauseTime - If the pause key is pressed, the duration of the pause is
% recorded to make sure we account for it when starting again!
function [ key, Resp_Time ] = getInput()

global NO_KEY 

%%
key = NO_KEY;

[KeyIsDown, Resp_Time, Resp1] = KbCheck();

if KeyIsDown
    key = find(Resp1);
end

if key == NO_KEY
    Resp_Time = [];
end
end