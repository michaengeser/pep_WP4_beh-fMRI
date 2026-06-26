function d = get_object_sim(d, cfg)

if ~isfield(cfg, 'analysis_names'); cfg.analysis_names = {'typical'  'control'  'photos'};end
if ~isfield(cfg, 'force_recompute'); cfg.force_recompute = false;end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true;end

warning off
% Computes pairwise object overlap and object sparsity comparison across all rows in the object CSV.
%
% Expected CSV format:
% id, object_1, object_2, object_3, ...
%
% Each row is treated as a SET of objects.
% Repeated objects within the same row are counted only once for object overlap.
% loop through categories

for category = cfg.categories
    category = char(category);

    for analysis_name = cfg.analysis_names
        analysis_name = char(analysis_name);

        %% Choose input file
        if strcmp(analysis_name, 'Own_drawing_DNN_RDM')
            filename = fullfile(pwd, '..', 'image_similarities', 'drawings_draw3D', ...
                'used_drawings_exp1', 'own', category, [category, '_objects.csv']);
        elseif strcmp(analysis_name, 'Control_drawing_DNN_RDM')
            filename = fullfile(pwd, '..', 'image_similarities', 'drawings_draw3D', ...
                'used_drawings_exp1', 'control', category, [category, '_control_objects.csv']);
        elseif strcmp(analysis_name, 'Photo_DNN_RDM')
            filename = fullfile(pwd, '..', 'image_similarities', 'pictures', category, [category, '_objects.csv']);
        end

        %% Read CSV
        T = readtable(filename, 'TextType', 'string', 'Delimiter', ',');

        % Filter for current subjects
        try
            T = T(ismember(T.id, cfg.subNums), :);
        catch
            [~,ids] = fileparts(T.id);
            for iid = 1:length(ids)
                idi = char(ids(iid));
                ids(iid) = idi(1:3);
            end
            T.id = str2double(ids);
            T = T(ismember(T.id, cfg.subNums), :);
        end

        % combine rows of the same subject
        if height(T) > cfg.n
            unique_ids = unique(T.id);
            T_all = T;
            T = T(1:cfg.n, :);
            T(:, 2:end) = []; % overwrite old table

            for iid = 1:cfg.n
                T.id(iid) = unique_ids(iid);
                sub_objs = T_all(T_all.id == unique_ids(iid),...
                    2:end);
                sub_objs = unique(table2array(sub_objs));
                for sub_obj = 1:length(sub_objs)
                    T{iid, 1 + sub_obj} = sub_objs(sub_obj);
                end 
            end 

        end 

        % First column should contain row IDs
        rowIDs = T{:, 1};

        % Remaining columns contain object names
        objectData = T{:, 2:end};

        nRows = height(T);

        %% Get unique object vocabulary across the whole file
        allObjects = objectData(:);

        % Remove missing/empty cells
        allObjects = strtrim(allObjects);
        allObjects = allObjects(~ismissing(allObjects) & allObjects ~= "");

        uniqueObjects = unique(allObjects);
        nObjects = numel(uniqueObjects);

        %% Build binary presence/absence matrix
        % rows = scenes/items
        % columns = unique objects
        presence = false(nRows, nObjects);
        subObjectCount = nan(1, cfg.n);

        for i = 1:nRows
            rowObjects = objectData(i, :);
            rowObjects = strtrim(rowObjects);
            rowObjects = rowObjects(~ismissing(rowObjects) & rowObjects ~= "");

            % Treat as a set: remove duplicates within row
            rowObjects = unique(rowObjects);
            subObjectCount(i) = numel(rowObjects);

            [tf, loc] = ismember(rowObjects, uniqueObjects);
            presence(i, loc(tf)) = true;
        end

        % get inter-subject RDM for object count
        objectCount_mat = abs(subObjectCount - subObjectCount');

        %% Compute pairwise object overlap and object sparsity
        objectSim_mat = zeros(nRows, nRows);

        for i = 1:nRows
            for j = 1:nRows
                intersectionCount = sum(presence(i, :) & presence(j, :));
                unionCount = sum(presence(i, :) | presence(j, :));

                if unionCount == 0
                    objectSim_mat(i, j) = NaN;
                else
                    objectSim_mat(i, j) = intersectionCount / unionCount;
                end
            end
        end

        if cfg.dissimilarity
            objectSim_mat = 1 - objectSim_mat;
        end

        objectSim.(analysis_name).(category).subject_mean.name = 'objectSim';
        objectSim.(analysis_name).(category).subject_mean.color = [0,0,0];
        objectSim.(analysis_name).(category).subject_mean.RDM = objectSim_mat;
        objectSim.(analysis_name).(category).subjects = rowIDs;

        objectCount.(analysis_name).(category).subject_mean.name = 'objectCount';
        objectCount.(analysis_name).(category).subject_mean.color = [0,0,0];
        objectCount.(analysis_name).(category).subject_mean.RDM = objectCount_mat;
        objectCount.(analysis_name).(category).subjects = rowIDs;

        disp("Computed pairwise object overlap similarity matrix for:");
        disp([category, ' - ', analysis_name]);

    end
end

% Write to results struct
d.DNN.objectSim = objectSim; % pretend this is a DNN approach
d.DNN.objectCount = objectCount; % pretend this is a DNN approach

disp("Done.");
warning on

end