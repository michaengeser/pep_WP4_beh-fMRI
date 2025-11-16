% need to be done

function [ ] = instructions(task)

global sub_num
global RightKey LeftKey

instr_slide1 = "instructions_slide1.png";

% decide which slides to take (based on task and subject number)
if strcmp(task, 'categorization')
    if mod(sub_num, 2)
        % odd subject number
        instr_slide2 = "instructions_cate_odd.png";
    else
        % even subject number
        instr_slide2 = "instructions_cate_even.png";
    end

    % last slide
    instr_slide3 = "instructions_end_cate.png";
else
    if strcmp(task, 'typicality')
        instr_slide2 = "instructions_typ.png";
    elseif strcmp(task, 'familiarity')
        instr_slide2 = "instructions_fam.png";
    elseif strcmp(task, 'usability')
        instr_slide2 = "instructions_usa.png";
    elseif strcmp(task, 'aesthetic')
        instr_slide2 = "instructions_aes.png";
    end

    % last slide
    instr_slide3 = "instructions_end_rating.png";
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
