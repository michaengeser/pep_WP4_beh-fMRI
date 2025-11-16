%% Housekeeping:
% Clearing the command window before we start saving it
sca;
close all;
clear all;

% To get different seeds for matlab randomization functions.
rng('shuffle');

% global parameters:
global sub_num TRUE FALSE refRate task SHOW_PRACTICE  session
global FRAME_ANTICIPATION PHOTODIODE DIOD_DURATION SHOW_INSTRUCTIONS category
global RESTART_KEY NO_KEY abortKey spaceBar valid_resp_keys
global expDir

expDir = pwd;

% Eye-tracking parameters
% global el EYE_TRACKER CalibrationKey ValidationKey EYETRACKER_CALIBRATION_MESSAGE compKbDevice

% Add functions folder to path (when we separate all functions)
function_folder = [pwd,filesep,'functions\'];
addpath(function_folder)

%% User input

% prompt user for information
sub_num = input('Subject number [101-199, default: 101]: ');
if isempty(sub_num); sub_num = 101; end % default
session = input('Session number [default: 1]: ');
if isempty(session); session = 1; end % default
task_num = input(['Task [0: Slide show, 1: Categorization, 2: Typicality, 3: Familiarity,', ...
    ' 4: Aesthetic, 5: Usability, 6: Complexity, default: 1]: ']);
if isempty(task_num); task_num = 1; end % default

% tanslate used input
task_array = {'slide_show', 'categorization', 'typicality', 'familiarity',...
    'aesthetic', 'usability', 'complexity'};
task = task_array{task_num + 1};

% ask for category if its the rating task
if task_num > 1
    category_num = input('Category [0: Bathroom, 1: Kitchen, default: 1]: ');
    if isempty(category_num); category_num = 1; end % default

    % Set the catagories:
    categories = ["bathroom", "kitchen"];

    % translate user input
    category = categories{category_num + 1};

else 
    category = ' ';
end

% Logging everything that is printed into the command window
dfile ='log_wp4_beh.txt';
if exist(dfile, 'file') ; delete(dfile); end
Str = CmdWinTool('getText');
dlmwrite(dfile,Str,'delimiter','');


%% check if participant and session exists already

SubSesFolder = fullfile(pwd,'data',['sub-', num2str(sub_num)],['ses-',num2str(session)], task);
ExistFlag = exist(SubSesFolder,'dir');
if ExistFlag
    WaitSecs(1);
    warning ('This participant number and session was already attributed!')
    proceedInput = questdlg({'This participant number and session was already attributed!', 'Are you sure you want to proceed?'},'RestartPrompt','yes','no','yes');
    if strcmp(proceedInput,'no')
        error('Program aborted by user')
    end
end

%% Initializing

% experimental parameters
initRuntimeParameters

% PTB:
initPsychtooblox(); % initializes psychtoolbox window at correct resolution and refresh rate

%% Load and prepare stimuli:
showMessage('Loading...');

%% Setup the trial matrix and log:
[tr_mat, file_list] = load_trial_matrix(sub_num);

% get trial matrix of the task
if strcmp(task, 'categorization')
    task_mat = tr_mat(strcmp(tr_mat.task, task),:);
else % if rating task, just take corresponding category
    task_mat = tr_mat(strcmp(tr_mat.task, task) & strcmp(tr_mat.category, category),:);
end

% get stimuli
loadStimuli(file_list)

% get mask
mask_texture = loadMask();

% make durations multiple of refresh rate
[task_mat] = durationXrefRate(task_mat);

%% Instructions
% displays instructions
if SHOW_INSTRUCTIONS
    instructions(task);
    WaitSecs(1);
end

%% Main experimental loop:
try

    %% save everything from command window
    Str = CmdWinTool('getText');
    dlmwrite(dfile,Str,'delimiter','');

    %% Block loop:
    blks = unique(task_mat.block);
    showFixation('PhotodiodeOff');

    if SHOW_PRACTICE
        blk = 0;
    else
        blk = 1;
    end

    while blk <= blks(end)

        %
        %         % Initialize the eyetracker with the block number and run the
        %         % calibration:
        %
        %         if EYE_TRACKER
        %             % Initialize the eyetracker:
        %             initEyetracker(sub_num, blk);
        %             % Show the calibration message to give the option to perform
        %             % the eyetracker calibration if needed:
        %             showMessage(EYETRACKER_CALIBRATION_MESSAGE);
        %             CorrectKey = 0; % Setting the CorrectKey to 0 to initiate the loop
        %             while ~CorrectKey % As long as a non-accepted key is pressed, keep on asking
        %                 [~, CalibrationResp, ~] = KbWait(compKbDevice,3);
        %                 if CalibrationResp(CalibrationKey)
        %                     % Run the calibration:
        %                     EyelinkDoTrackerSetup(el);
        %                     CorrectKey = 1;
        %                 elseif CalibrationResp(ValidationKey)
        %                     CorrectKey = 1;
        %                 end
        %             end
        %             % Starting the recording
        %             Eyelink('StartRecording');
        %             % Wait for the recording to have started:
        %             WaitSecs(0.1);
        %         end

        % Extract the trial and log of this block only:
        blk_mat = task_mat(task_mat.block == blk, :);

        % Add the columns for logging:
        blk_mat = prepare_log(blk_mat);

        % check what block we´re in and informs user
        if blk == 0
            showMessage('Press a space to start the practice');
        elseif blk == 1
            showMessage('Press a space to start the experiment');
        else
            showMessage('Press a space to start the next block');
        end

        % wait for space bar
        [~, ~, wait_resp] = KbCheck();
        while ~wait_resp(spaceBar)
            [~, ~, wait_resp] = KbCheck();
        end

        % Wait a random amount of time and show fixation:
        showFixation('PhotodiodeOff');
        WaitSecs(rand + 2);

        %% Trials loop:
        for tr = 1:length(blk_mat.trial)

            % flags needs to be initialized
            fixShown = FALSE;
            maskShown = FALSE;
            jitterLogged = FALSE;
            hasInput = FALSE;

            % show stimulus
            blk_mat.stim_time(tr) = showStimuli(blk_mat.texture(tr));
            DiodFrame = 0;

            %             % Sending response trigger for the eyetracker
            %             if EYE_TRACKER
            %                 trigger_str = get_et_trigger('vis_onset', blk_mat.task_relevance{tr}, ...
            %                     blk_mat.duration(tr), blk_mat.category{tr}, orientation, vis_stim_id, ...
            %                     blk_mat.SOA(tr), blk_mat.SOA_lock(tr), blk_mat.pitch(tr));
            %                 Eyelink('Message',trigger_str);
            %             end



            %--------------------------------------------------------

            %% TIME LOOP

            % I then set a frame counter. The flip of the stimulus
            % presentation is frame 0. It is already the previous frame because it already occured:
            PreviousFrame = 0;
            % I then set a frame index. It is the same as the previous
            % frame for now
            FrameIndex = PreviousFrame;

            % start timer
            elapsedTime = 0;

            % define total trial duration
            min_trial_duration = blk_mat.duration(tr) + blk_mat.mask_dur(tr);

            while elapsedTime < min_trial_duration || ~hasInput

                %% Get response:
                if ~hasInput
                    % Ge the response:
                    [key,Resp_Time] = getInput();

                    % Handling the response:
                    if key ~= NO_KEY

                        %                         % Sending response trigger for the eyetracker
                        %                         if EYE_TRACKER
                        %                             trigger_str = get_et_trigger('response', blk_mat.task_relevance{tr}, ...
                        %                                 blk_mat.duration(tr), blk_mat.category{tr}, orientation, vis_stim_id, ...
                        %                                 blk_mat.SOA(tr), blk_mat.SOA_lock(tr), blk_mat.pitch(tr));
                        %                             Eyelink('Message',trigger_str);
                        %                         end

                        % If the experiment was aborted:
                        if key == abortKey
                            error('Experiment has been aborted');

                            % if pressed key is not a valid response key
                        elseif ~ismember(key, valid_resp_keys)
                            invalid_key_msg = ['Invalid response key!', newline, newline,...
                                'Press space to proceed'];
                            showMessage(invalid_key_msg);

                            % wait for space bar
                            [~, ~, wait_resp] = KbCheck();
                            while ~wait_resp(spaceBar)
                                [~, ~, wait_resp] = KbCheck();
                            end

                        end

                        % Log the response received:
                        blk_mat.trial_response(tr) = key;
                        blk_mat.time_of_resp(tr) =  Resp_Time;

                        % logging reaction
                        hasInput = TRUE;
                    end
                end

                %% mask

                % check if trial has a mask
                if blk_mat.mask_dur(tr) > 0

                    % Present mask
                    if elapsedTime >= (blk_mat.duration(tr) - refRate*FRAME_ANTICIPATION) && maskShown == FALSE
                        mask_time = showStimuli(mask_texture);
                        DiodFrame = CurrentFrame;

                        %                     % Sending response trigger for the eyetracker
                        %                     if EYE_TRACKER
                        %                         trigger_str = get_et_trigger('mask_onset', blk_mat.task_relevance{tr}, ...
                        %                             blk_mat.duration(tr), blk_mat.category{tr}, orientation, vis_stim_id, ...
                        %                             blk_mat.SOA(tr), blk_mat.SOA_lock(tr), blk_mat.pitch(tr));
                        %                         Eyelink('Message',trigger_str);
                        %                     end

                        % log fixation
                        blk_mat.mask_time(tr) = mask_time;
                        maskShown = TRUE;
                    end
                end

                %% Waiting for response

                % Present fixation while waiting for response
                if elapsedTime >= (min_trial_duration - refRate*FRAME_ANTICIPATION) && fixShown == FALSE
                    fix_time = showFixation('PhotodiodeOn');
                    DiodFrame = CurrentFrame;

                    %                     % Sending response trigger for the eyetracker
                    %                     if EYE_TRACKER
                    %                         trigger_str = get_et_trigger('fixation_onset', blk_mat.task_relevance{tr}, ...
                    %                             blk_mat.duration(tr), blk_mat.category{tr}, orientation, vis_stim_id, ...
                    %                             blk_mat.SOA(tr), blk_mat.SOA_lock(tr), blk_mat.pitch(tr));
                    %                         Eyelink('Message',trigger_str);
                    %                     end

                    % log fixation
                    blk_mat.fix_time(tr) = fix_time;
                    fixShown = TRUE;
                end

                % Updating clock:
                elapsedTime = GetSecs - blk_mat.stim_time(tr);

                % Updating the frame counter:
                CurrentFrame = floor(elapsedTime/refRate);

                % Check if a new frame started:
                if CurrentFrame > PreviousFrame
                    FrameIndex = FrameIndex +1;

                    % turn photodiode off again after diod duration
                    if PHOTODIODE && (CurrentFrame - DiodFrame == DIOD_DURATION - 1)
                        turnPhotoTrigger('off');
                    end
                    PreviousFrame = CurrentFrame;
                end
            end

            %% Jitter TIME LOOP

            % start timer
            elapsedTime = 0;

            % define total trial duration
            jitter_duration = blk_mat.jitter(tr) - (refRate*FRAME_ANTICIPATION) ;

            while elapsedTime < jitter_duration

                % Present jitter
                if jitterLogged == FALSE

                    blk_mat.jit_onset(tr) = showFixation('PhotodiodeOn');
                    DiodFrame = CurrentFrame;

                    %                     % Sending response trigger for the eyetracker
                    %                     if EYE_TRACKER
                    %                         trigger_str = get_et_trigger('jitter_onset', blk_mat.task_relevance{tr}, ...
                    %                             blk_mat.duration(tr), blk_mat.category{tr}, orientation, vis_stim_id, ...
                    %                             blk_mat.SOA(tr), blk_mat.SOA_lock(tr), blk_mat.pitch(tr));
                    %                         Eyelink('Message',trigger_str);
                    %                     end

                    jitterLogged = TRUE;

                end
                % Updating clock:
                elapsedTime = GetSecs - blk_mat.jit_onset(tr);

            end

            % log end of trial
            blk_mat.trial_end(tr) = GetSecs;

            %% End of trial

            if(key==RESTART_KEY)
                break
            end
        end  % End of trial loop

        % Save the data of this block:
        saveTable(blk_mat, blk);

        %         % Save the eyetracker data:
        %         if EYE_TRACKER
        %             saveEyetracker(task, blk);
        %         end

        % Append the block log to the overall log:
        if ~exist('log_all', 'var')
            log_all = blk_mat;
        else
            log_all = [log_all; blk_mat];  % Not the most efficient but it is in a non critical part
        end

        %% Feedback
        blk = blockFeedback(blk_mat, task_mat);

        % Move to the next block
        blk = blk + 1;
        WaitSecs(0.2);

    end  % End of block loop

    %% End of experiment

    % Letting the participant know that it is over:
    end_message = ['THE END', newline, newline,'Thank you!'];
    showMessage(end_message);

    %% Save

    % Save the whole table:
    saveTable(log_all, "all");

    % Save the code:
    saveCode(task);

    WaitSecs(2);
    showMessage('saving...');

    % save everything from command window
    Str = CmdWinTool('getText');
    dlmwrite(dfile,Str,'delimiter','');

    % Terminating teh experiment:
    safeExit()

catch e

    % Save the data:
    try
        % Save the beh data:
        saveTable(blk_mat, blk);

        % Save the eyetracker data:
        %         if EYE_TRACKER
        %             saveEyetracker(task, blk);
        %         end

        % If the log all already exists, save it as well:
        if exist('log_all', 'var')
            [log_all] = compute_performance(log_all);
            saveTable(log_all, "all");
        end

        % Save the code:
        saveCode(task);
        safeExit()

    catch ee
        warning('-----  Data could not be saved!  ------')
        safeExit()
        rethrow(ee);
    end

    rethrow(e);
end