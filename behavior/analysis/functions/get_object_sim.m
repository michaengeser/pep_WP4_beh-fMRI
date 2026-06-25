function d = get_object_sim(d, cfg)

if ~isfield(cfg, 'analysis_names'); cfg.analysis_names = {'typical'  'control'  'photos'};end
if ~isfield(cfg, 'force_recompute'); cfg.force_recompute = false;end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true;end

% Computes pairwise object overlap and object sparsity comparison across all rows in the object CSV.
%
% Expected CSV format:
% id, object_1, object_2, object_3, ...
%
% Each row is treated as a SET of objects.
% Repeated objects within the same row are counted only once for object overlap.
% loop through categories

outputFilename = fullfile(pwd, '..', 'derivatives', 'group_level', 'object_overlap.mat');

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
                    'used_drawings_exp1', 'own', category, [category, '_objects.csv']);
            elseif strcmp(analysis_name, 'control')
                filename = fullfile(pwd, '..', 'image_similarities', 'drawings_draw3D', ...
                    'used_drawings_exp1', 'control', category, [category, '_control_objects.csv']);
            elseif strcmp(analysis_name, 'photos')
                filename = fullfile(pwd, '..', 'image_similarities', 'pictures', category, [category, '_objects.csv']);
            end

            %% Read CSV
            T = readtable(filename, 'TextType', 'string');

            % Filter for current subjects
            T = T(ismember(T.id, cfg.subNums), :);

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

    % Save results
    save(outputFilename, "objectSim", "objectCount")
end

% Write to results struct
d.DNN.objectSim = objectSim; % pretend this is a DNN approach
d.DNN.objectCount = objectCount; % pretend this is a DNN approach

disp("Done.");

end