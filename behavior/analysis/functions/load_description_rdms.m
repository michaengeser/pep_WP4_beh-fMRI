function d = load_description_rdms(d, cfg)

warning('off')

if ~isfield(cfg, 'analysis_names'); cfg.analysis_names = {'typical'  'control'  'photos'};end
if ~isfield(cfg, 'force_recompute'); cfg.force_recompute = true;end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true;end

for category = cfg.categories
    category = char(category);

    for analysis_name = cfg.analysis_names
        analysis_name = char(analysis_name);

        %% Choose input file
        if strcmp(analysis_name, 'Own_drawing_DNN_RDM')
            filename = fullfile(pwd, '..', 'image_similarities', 'drawings_draw3D', ...
                'used_drawings_exp1', 'own', category,...
                ['gem_2-5_descriptions_', category, '_MPNet_spearman.csv']);
        elseif strcmp(analysis_name, 'Control_drawing_DNN_RDM')
            filename = fullfile(pwd, '..', 'image_similarities', 'drawings_draw3D', ...
                'used_drawings_exp1', 'control', category,...
                ['gem_2-5_descriptions_', category, '_control_MPNet_spearman.csv']);
        elseif strcmp(analysis_name, 'Photo_DNN_RDM')
            filename = fullfile(pwd, '..', 'image_similarities', 'pictures', category,...
                ['gem_2-5_descriptions_', category, '_photos_MPNet_spearman.csv']);
        end

        %% Read CSV
        T = readtable(filename, 'TextType', 'string');

        % Filter for current subjects
        is_rdm = zeros(cfg.n, cfg.n);
        row_names = T(:, 1);
        col_names = T.Properties.VariableNames;

        for iSub1 = 1:cfg.n
            sub1 = cfg.subNums(iSub1);
            row_idx = contains(table2array(row_names), num2str(sub1));

            for iSub2 = 1:cfg.n
                sub2 = cfg.subNums(iSub2);
                col_idx = contains(col_names, num2str(sub2));


                % get cells from subject
                is_rdm(iSub1, iSub2) = mean(table2array(T(row_idx, col_idx)), "all");
            end
        end

        descSim.(analysis_name).(category).subject_mean.name = 'descSim';
        descSim.(analysis_name).(category).subject_mean.color = [0,0,0];
        descSim.(analysis_name).(category).subject_mean.RDM = is_rdm;
        descSim.(analysis_name).(category).subjects = cfg.subNums;

        disp("Computed pairwise description similarity matrix for:");
        disp([category, ' - ', analysis_name]);

    end
end

% Write to results struct
d.DNN.descSim = descSim; % pretend this is a DNN approach

disp("Done.");
warning('on')
end