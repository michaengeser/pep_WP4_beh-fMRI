% need to be done

function [ ] = instructions()

global category task
global RightKey LeftKey

% decide which slides to take (based on task and subject number)
if strcmp(task, 'categorization')

    instr_slide1 = 'Slide01.tif';
    instr_slide2 = 'Slide02.tif';
    instr_slide3 = 'Slide03.tif';

else

    if strcmp(category, 'kitchen')
        instr_slide1 = 'Slide04.tif';
    elseif strcmp(category, 'bathroom')
        instr_slide1 = 'Slide05.tif';
    end

    if strcmp(task, 'typicality')
        instr_slide2 = 'Slide06.tif';
    elseif strcmp(task, 'familiarity')
        instr_slide2 = 'Slide07.tif';
    elseif strcmp(task, 'aesthetic')
        instr_slide2 = 'Slide08.tif';
    elseif strcmp(task, 'usability')
        instr_slide2 = 'Slide09.tif';
    elseif strcmp(task, 'complexity')
        instr_slide2 = 'Slide10.tif';
    end

    % last slide
    instr_slide3 = 'Slide11.tif';
end
        
instrsPaths = [instr_slide1; instr_slide2; instr_slide3];

% Setting the slide number to 1 to initiate the while loop
instrSlideNum = 1;

% displays all instr screens
while instrSlideNum<= size(instrsPaths,1) % Looping until we went through all slides:

    % Showing instr slide
    showInstructions(instrsPaths(instrSlideNum,:));
    WaitSecs(0.3);

    CorrectKey = 0; % Setting the CorrectKey to 0 to initiate the loop
    while ~CorrectKey % As long as a non-accepted key is pressed, keep on asking
        [~, ~, instrsResp] = KbCheck();

        if instrsResp(RightKey) % If the participant press the right key, increment by 1 the slide number
            instrSlideNum = instrSlideNum + 1;
            CorrectKey = 1;

        elseif instrsResp(LeftKey) % Else if the participant pressed the left key:
            if instrSlideNum == 1 % If we are at slide one, that doesn't work
                CorrectKey = 0;
            else % Otherwise, just go back one slide
                instrSlideNum = instrSlideNum - 1;
                CorrectKey = 1;
            end
        end
    end
end
end
