function d = get_image_betas(d, cfg)

if ~isfield(cfg, 'force_recompute'); cfg.force_recompute = false; end

sourcedataPath = fullfile(pwd, '..','sourcedata');

% check if file exists
file_name = fullfile(pwd, '..','derivatives', 'group_level', 'all_betas.mat');
if exist(file_name, 'file') && ~cfg.force_recompute
    disp('Loading existing data')
    load(file_name)
else

    %% make multiple condition files
    cfg.sortRows = true;
    cfg.includeTargets = true;
    create_mcf_func(cfg)

    %% load ROI mask
    nmasks = numel(cfg.rois);

    for iSub = 1:length(cfg.subNums)

        % subject ID
        subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
        subID2 = sprintf('sub%0.3d', cfg.subNums(iSub));

        % get image name and trial ID correspondence
        template_names = [];
        template_trialIDs = [];
        for i = 1:2
            mcf = fullfile(sourcedataPath, subID, 'beh', 'onsets', ...
                sprintf('mcf_%s_run-%s.mat', subID, num2str(i)));
            mcf_file = load(mcf);

            template_names = [template_names; unique(mcf_file.names)];
            template_trialIDs = [template_trialIDs; unique(cell2mat(mcf_file.trialIDs))];
        end

        % load GLM single betas
        GLMsinglePath = fullfile(pwd, '..', 'derivatives', subID, 'GLMsingleEstimates');
        betaPath = fullfile(GLMsinglePath, 'GLMsingle_betas.nii');

        % if betas do not exist, compute them
        if ~exist(betaPath, 'file')
            cfg_temp = cfg;
            cfg_temp.subNums = cfg.subNums(iSub);
            GLMsingle_new_script(cfg_temp)
        end

        % load beta map
        beta_img = load_untouch_nii(betaPath);
        beta_img = beta_img.img;

        % get trial IDs after GLM single
        load(fullfile(GLMsinglePath, 'trialIDs.mat'));

        %% average across voxels for each mask
        for j=1:nmasks

            % get mask name
            mask_label=cfg.rois{j};
            mask_label_short = split(mask_label, '.');
            mask_label_short = mask_label_short{1};
            mask_label_short = mask_label_short(2:end);

            % progress report
            disp(['Extracting betas for: ', subID, ' - ', mask_label_short])


            % check if functional or anatomical ROI
            if ismember(mask_label_short, {'PPA', 'TOS', 'RSC', 'LOC', 'LPFC'})
                mask_fn=fullfile(pwd, '..', 'MNI_ROIs', 'func_ROIs', subID,...
                    [mask_label_short, '_funcROI.nii']);
            else
                mask_fn=fullfile(pwd, '..', 'MNI_ROIs', [char(mask_label)]);
            end

            % get mask
            mask = load_untouch_nii(mask_fn);
            newMaskImg =double(mask.img);
            if max(max(max(double(newMaskImg)))) > 1
                newMaskImg = newMaskImg/max(max(max(newMaskImg)));
            end

            % make 4D mask
            currentROI = newMaskImg;
            currentROI4D = repmat(currentROI, ...
                [1, 1, 1, size(beta_img, 4)]);

            % filter for data and take mean
            currentImage = beta_img;
            filtered = currentImage(currentROI4D == 1);
            ROIImage = reshape(filtered, [], size(beta_img, 4));
            meanData = mean(ROIImage, 1, 'omitnan');

            % loop trial IDs
            for trID = template_trialIDs'
                if trID <= cfg.nTrials/2
                    category = cfg.categories{1};
                    idx = trID;
                else
                    category = cfg.categories{2};
                    idx = trID - cfg.nTrials/2;
                end

                betas.(category).(subID2).(mask_label_short)(idx).image_name = template_names{template_trialIDs == trID};
                betas.(category).(subID2).(mask_label_short)(idx).all_betas = meanData(trID == trialIDs(:, 1));
                betas.(category).(subID2).(mask_label_short) (idx).mean_betas = mean(meanData(trID == trialIDs(:, 1)));
            end
        end % rois
    end

    % save 
    save(file_name, 'betas')
end

% store results in d strucure 
d.betas = betas;

end