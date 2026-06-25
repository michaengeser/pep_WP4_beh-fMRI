function d = compute_spatial_coherence_rdms(d, cfg)

warning('off')

if ~isfield(cfg, 'analysis_names'); cfg.analysis_names = {'typical'  'control'  'photos'};end
if ~isfield(cfg, 'force_recompute'); cfg.force_recompute = false;end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true;end

outputFilename = fullfile(pwd, '..', 'derivatives', 'group_level', 'sc_ce_similarity.mat');

if exist(outputFilename, 'file') && cfg.force_recompute == false
    load(outputFilename)
else

    for category = cfg.categories
        category = char(category);

        for analysis_name = cfg.analysis_names
            analysis_name = char(analysis_name);

            %% Choose input file
            if strcmp(analysis_name, 'typical')
                filename = fullfile(pwd, '..', 'image_similarities', 'drawings_draw3D', ...
                    'draw3D_images_exp1', 'own', category, 'LGNstatistics.mat');
                gen_tag = 'gen';
                copy_tag = 't_g'; % copy files do not have this (i.e. 102_bat_copy_gen3.png)
            elseif strcmp(analysis_name, 'control')
                filename = fullfile(pwd, '..', 'image_similarities', 'drawings_draw3D', ...
                    'draw3D_images_exp1', 'control', category, 'LGNstatistics.mat');
                gen_tag = 'gen';
                copy_tag = 'copy';
            elseif strcmp(analysis_name, 'photos')
                filename = fullfile(pwd, '..', 'image_similarities', 'pictures', category, 'LGNstatistics.mat');
                gen_tag = '';
                copy_tag = '';
            end

            %% Read data
            load(filename);

            % Filter for current images
            current_idx = (contains(filenames, category(1:3)) & ...
                contains(filenames, gen_tag) & ...
                contains(filenames, copy_tag));
            current_names = filenames(current_idx);
            current_SC = SC(current_idx, 1);
            current_CE = CE(current_idx, 1);

            % get inter-subject RDM
            sub_sc = zeros(1, cfg.n);
            sub_ce = zeros(1, cfg.n);
            for iSub = 1:cfg.n
                sub_sc(iSub) = mean(current_SC(contains(current_names, num2str(cfg.subNums(iSub)))));
                sub_ce(iSub) = mean(current_CE(contains(current_names, num2str(cfg.subNums(iSub)))));
            end 
            sc_rdm = abs(sub_sc - sub_sc');
            ce_rdm = abs(sub_ce - sub_ce');

            sc.(analysis_name).(category).subject_mean.name = 'sc';
            sc.(analysis_name).(category).subject_mean.color = [0,0,0];
            sc.(analysis_name).(category).subject_mean.RDM = sc_rdm;
            sc.(analysis_name).(category).subjects = cfg.subNums;

            ce.(analysis_name).(category).subject_mean.name = 'ce';
            ce.(analysis_name).(category).subject_mean.color = [0,0,0];
            ce.(analysis_name).(category).subject_mean.RDM = ce_rdm;
            ce.(analysis_name).(category).subjects = cfg.subNums;

            disp("Computed pairwise contrast similarity matrix for:");
            disp([category, ' - ', analysis_name]);

        end
    end

    % Save results
    save(outputFilename, "sc", "ce")
end

% Write to results struct
d.DNN.sc = sc; % pretend this is a DNN approach
d.DNN.ce = ce; % pretend this is a DNN approach

disp("Done.");
warning('on')
end