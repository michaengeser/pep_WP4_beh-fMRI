% this function makes a quick analysisfor a single subject of the
% reconstructed time prp task. As inout it recieves the subject number, the
% sessions and the task type. It returns a structure that contains the
% performances and some controls for the presentation of the stimuli

function [output_struct,refRate] = quick_analysis(sub_num, session, task, plotting)

if nargin<2
    session = 1;
    task = 'categorization';
    plotting = 0;
elseif nargin<3
    task = 'categorization';
    plotting = 0;
elseif nargin<4
    plotting = 0;
end
%%
global refRate
if isempty(refRate); refRate = 1/60; end

% make file name
filename = string(fullfile(pwd,'data',['sub-', num2str(sub_num)],['ses-',num2str(session)],task,...
    ['sub-', num2str(sub_num),'_ses-',num2str(session),'_run-all_task-', task,'_events.mat']));

event_table = load(filename);

try
    event_table = event_table.input_table;
catch
    event_table = event_table.d;
end

if ~ ismember('trial_accuracy', event_table.Properties.VariableNames)
    % Add functions folder to path (when we separate all functions)
    function_folder = [parent_dir,filesep,'functions\'];
    addpath(function_folder)

    [event_table, ~] = compute_performance(event_table);

end

% remove targets and practice
practice_event_table = event_table(event_table.is_practice == 1, :);
full_event_table = event_table;
% event_table = event_table(~event_table.is_practice,:);


% %% reaction time
% 
% % reaction time
% output_struct.mean_RT_aud = mean(event_table.RT_aud, 'omitnan');
% 
% % split by SOA and SOA lock
% output_struct.RT_aud_off(1) = mean(event_table.RT_aud(event_table.SOA == 0 & strcmp(event_table.SOA_lock, 'offset')), 'omitnan');
% output_struct.RT_aud_off(2) = mean(event_table.RT_aud(event_table.SOA == 0.116 & strcmp(event_table.SOA_lock, 'offset')), 'omitnan');
% output_struct.RT_aud_off(3) = mean(event_table.RT_aud(event_table.SOA == 0.232 & strcmp(event_table.SOA_lock, 'offset')), 'omitnan');
% output_struct.RT_aud_off(4) = mean(event_table.RT_aud(event_table.SOA == 0.466 & strcmp(event_table.SOA_lock, 'offset')), 'omitnan');
% 
% output_struct.RT_aud_on(1) = mean(event_table.RT_aud(event_table.SOA == 0 & strcmp(event_table.SOA_lock, 'onset')), 'omitnan');
% output_struct.RT_aud_on(2) = mean(event_table.RT_aud(event_table.SOA == 0.116 & strcmp(event_table.SOA_lock, 'onset')), 'omitnan');
% output_struct.RT_aud_on(3) = mean(event_table.RT_aud(event_table.SOA == 0.232 & strcmp(event_table.SOA_lock, 'onset')), 'omitnan');
% output_struct.RT_aud_on(4) = mean(event_table.RT_aud(event_table.SOA == 0.466 & strcmp(event_table.SOA_lock, 'onset')), 'omitnan');
% 
% % split by task relevnace
% 
% % task relevant
% relevant_event_table = event_table(strcmp(event_table.task_relevance, 'non-target'),:);
% 
% output_struct.task_relevant.RT_aud_off(1) = mean(relevant_event_table.RT_aud(relevant_event_table.SOA == 0 & strcmp(relevant_event_table.SOA_lock, 'offset')), 'omitnan');
% output_struct.task_relevant.RT_aud_off(2) = mean(relevant_event_table.RT_aud(relevant_event_table.SOA == 0.116 & strcmp(relevant_event_table.SOA_lock, 'offset')), 'omitnan');
% output_struct.task_relevant.RT_aud_off(3) = mean(relevant_event_table.RT_aud(relevant_event_table.SOA == 0.232 & strcmp(relevant_event_table.SOA_lock, 'offset')), 'omitnan');
% output_struct.task_relevant.RT_aud_off(4) = mean(relevant_event_table.RT_aud(relevant_event_table.SOA == 0.466 & strcmp(relevant_event_table.SOA_lock, 'offset')), 'omitnan');
% 
% output_struct.task_relevant.RT_aud_on(1) = mean(relevant_event_table.RT_aud(relevant_event_table.SOA == 0 & strcmp(relevant_event_table.SOA_lock, 'onset')), 'omitnan');
% output_struct.task_relevant.RT_aud_on(2) = mean(relevant_event_table.RT_aud(relevant_event_table.SOA == 0.116 & strcmp(relevant_event_table.SOA_lock, 'onset')), 'omitnan');
% output_struct.task_relevant.RT_aud_on(3) = mean(relevant_event_table.RT_aud(relevant_event_table.SOA == 0.232 & strcmp(relevant_event_table.SOA_lock, 'onset')), 'omitnan');
% output_struct.task_relevant.RT_aud_on(4) = mean(relevant_event_table.RT_aud(relevant_event_table.SOA == 0.466 & strcmp(relevant_event_table.SOA_lock, 'onset')), 'omitnan');
% 
% 
% % task irrelevant
% irrelevant_event_table = event_table(strcmp(event_table.task_relevance, 'irrelevant'),:);
% 
% output_struct.task_irrelevant.RT_aud_off(1) = mean(irrelevant_event_table.RT_aud(irrelevant_event_table.SOA == 0 & strcmp(irrelevant_event_table.SOA_lock, 'offset')), 'omitnan');
% output_struct.task_irrelevant.RT_aud_off(2) = mean(irrelevant_event_table.RT_aud(irrelevant_event_table.SOA == 0.116 & strcmp(irrelevant_event_table.SOA_lock, 'offset')), 'omitnan');
% output_struct.task_irrelevant.RT_aud_off(3) = mean(irrelevant_event_table.RT_aud(irrelevant_event_table.SOA == 0.232 & strcmp(irrelevant_event_table.SOA_lock, 'offset')), 'omitnan');
% output_struct.task_irrelevant.RT_aud_off(4) = mean(irrelevant_event_table.RT_aud(irrelevant_event_table.SOA == 0.466 & strcmp(irrelevant_event_table.SOA_lock, 'offset')), 'omitnan');
% 
% output_struct.task_irrelevant.RT_aud_on(1) = mean(irrelevant_event_table.RT_aud(irrelevant_event_table.SOA == 0 & strcmp(irrelevant_event_table.SOA_lock, 'onset')), 'omitnan');
% output_struct.task_irrelevant.RT_aud_on(2) = mean(irrelevant_event_table.RT_aud(irrelevant_event_table.SOA == 0.116 & strcmp(irrelevant_event_table.SOA_lock, 'onset')), 'omitnan');
% output_struct.task_irrelevant.RT_aud_on(3) = mean(irrelevant_event_table.RT_aud(irrelevant_event_table.SOA == 0.232 & strcmp(irrelevant_event_table.SOA_lock, 'onset')), 'omitnan');
% output_struct.task_irrelevant.RT_aud_on(4) = mean(irrelevant_event_table.RT_aud(irrelevant_event_table.SOA == 0.466 & strcmp(irrelevant_event_table.SOA_lock, 'onset')), 'omitnan');


%% controls

% control stim duration
output_struct.real_dur = event_table.mask_time - event_table.stim_time;
output_struct.dur_diff = output_struct.real_dur - event_table.duration;
output_struct.dur_diff_max = max(output_struct.dur_diff);
output_struct.dur_diff_mean = mean(output_struct.dur_diff);

% check currect frames
output_struct.correct_frame = sum(output_struct.dur_diff < (0.5*refRate) & output_struct.dur_diff > (-0.5*refRate));
output_struct.one_frame_late = sum(output_struct.dur_diff > (0.5*refRate) & output_struct.dur_diff < (1.5*refRate));
output_struct.two_frames_late = sum(output_struct.dur_diff > (1.5*refRate) & output_struct.dur_diff < (2.5*refRate));
output_struct.three_frames_late = sum(output_struct.dur_diff > (2.5*refRate) & output_struct.dur_diff < (3.5*refRate));
output_struct.one_frame_early = sum(output_struct.dur_diff < (-0.5*refRate) & output_struct.dur_diff > (-1.5*refRate));
output_struct.two_frames_early = sum(output_struct.dur_diff < (-1.5*refRate) & output_struct.dur_diff > (-2.5*refRate));
output_struct.three_frames_early = sum(output_struct.dur_diff < (-2.5*refRate) & output_struct.dur_diff > (-3.5*refRate));

% control trial duration
output_struct.trial_durs = event_table.trial_end - event_table.stim_time;
output_struct.trial_dur_mean = mean(output_struct.trial_durs);

output_struct.target_trial_durs = event_table.RT + event_table.jitter;
output_struct.diff_trial_durs = output_struct.trial_durs - output_struct.target_trial_durs;

% control masking + stimulus timing
output_struct.ctr_mask_real_dur = event_table.fix_time - event_table.stim_time;
output_struct.ctr_mask_dur_diff = output_struct.ctr_mask_real_dur - (event_table.duration + event_table.mask_dur);
output_struct.ctr_mask_dur_diff_max = max(output_struct.ctr_mask_dur_diff);
output_struct.ctr_mask_dur_diff_mean = mean(output_struct.ctr_mask_dur_diff);


% control masking timing
output_struct.mask_real_dur = event_table.fix_time - event_table.mask_time;
output_struct.mask_dur_diff = output_struct.mask_real_dur - (event_table.mask_dur);
output_struct.mask_dur_diff_max = max(output_struct.mask_dur_diff);
output_struct.mask_dur_diff_mean = mean(output_struct.mask_dur_diff);



%% timing control
% 
% % timing control table
% output_timing_table = table;
% output_timing_table.real_dur = output_struct.real_dur;
% output_timing_table.target_dur = event_table.duration;
% output_timing_table.dur_diff = output_struct.dur_diff;
% output_timing_table.real_SOA = output_struct.real_SOA;
% output_timing_table.target_SOA = event_table.onset_SOA;
% output_timing_table.SOA_diff = output_struct.SOA_diff;
% output_timing_table.trial_dur = output_struct.trial_durs;
% output_timing_table.real_trial_dur_plus_jitter = output_struct.real_trial_durs_plus_jitter;
% output_timing_table.taregt_trial_dur_plus_jitter = output_struct.target_trial_durs_plus_jitter;
% output_timing_table.diff_trial_dur_plus_jitter = output_struct.diff_trial_durs_plus_jitter;
% 
% timing_output_filename = string(fullfile(test_data_dir, ['sub-', num2str(sub_num),'_ses-',num2str(session),'_', task,'_test_timing.csv']));
% writetable(output_timing_table, timing_output_filename);

%% plotting

% if plotting
% 
%     %% timing plots
%     figure(2)
%     histogram(output_struct.dur_diff, [-3.5*refRate, -2.5*refRate, -1.5*refRate, -0.5*refRate, ...
%         0.5*refRate, 1.5*refRate, 2.5*refRate, 3.5*refRate]);
%     ylim([0, length(output_struct.dur_diff) + 50])
%     title(['Jitter in visual stimulus duration (Refresh rate: ', num2str(refRate), ' s)'])
%     ylabel('Number of trials')
%     xlabel('Frames off')
% 
%     h1 = gca;
%     h1.XTick = [-3*refRate, -2*refRate, -refRate, 0, refRate, 2*refRate, 3*refRate];
%     h1.XTickLabel = ["-3", "-2", "-1", "0", "+1", "+2", "+3"];
% 
%     fig_filename = string(fullfile(test_data_dir, ['sub-', num2str(sub_num),'_ses-',num2str(session),'_', task,'_duration_jitter.fig']));
%     saveas(h1,fig_filename)
%     png_filename = string(fullfile(test_data_dir, ['sub-', num2str(sub_num),'_ses-',num2str(session),'_', task,'_duration_jitter.png']));
%     saveas(h1,png_filename)
% 
%     figure(3)
%     histogram(output_struct.SOA_diff*1000);
%     title('Jitter in SOA')
%     ylabel('Number of trials')
%     xlabel('ms off')
% 
%     h2 = gca;
%     fig_filename = string(fullfile(test_data_dir, ['sub-', num2str(sub_num),'_ses-',num2str(session),'_', task,'_all_SOA_jitter.fig']));
%     saveas(h2,fig_filename)
%     png_filename = string(fullfile(test_data_dir, ['sub-', num2str(sub_num),'_ses-',num2str(session),'_', task,'_all_SOA_jitter.png']));
%     saveas(h2,png_filename)
% 
%     figure(4)
%     single_soas = tiledlayout(2,2);
%     nbins = 8;
%     for sub_plot = 1:4
%         nexttile
% 
%         onset_data = (output_struct.SOA_diff(event_table.SOA == SOAs(sub_plot) & strcmp(event_table.SOA_lock, 'onset')))*1000;
%         offset_data = (output_struct.SOA_diff(event_table.SOA == SOAs(sub_plot) & strcmp(event_table.SOA_lock, 'offset')))*1000;
%         plot_data = nan(length(offset_data)*2, 2);
%         plot_data(1:length(onset_data),1) = onset_data;
%         plot_data(1:length(offset_data),2) = offset_data;
%         hist(plot_data,nbins);
%         title(['SOA: ', num2str(SOAs(sub_plot)), ' ms'])
%         ylabel('Number of trials')
%         xlabel('ms off')
%     end
% 
%     lgd = legend('offset', 'onset');
%     lgd.Layout.Tile = 'east';
% 
%     fig_filename = string(fullfile(test_data_dir, ['sub-', num2str(sub_num),'_ses-',num2str(session),'_', task,'_single_SOAs_jitter.fig']));
%     saveas(single_soas,fig_filename)
%     png_filename = string(fullfile(test_data_dir, ['sub-', num2str(sub_num),'_ses-',num2str(session),'_', task,'_single_SOAs_jitter.png']));
%     saveas(single_soas,png_filename)
% 
% end
end