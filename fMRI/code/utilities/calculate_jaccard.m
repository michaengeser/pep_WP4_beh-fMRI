function d = calculate_jaccard(d, cfg)

if ~isfield(cfg, 'analysis_names'); cfg.analysis_names = {'typical'  'control'  'photos'};end
if ~isfield(cfg, 'force_recompute'); cfg.force_recompute = true;end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true;end

% calculate_jaccard_from_csv.m
% Computes pairwise Jaccard similarity across all rows in the object CSV.
%
% Expected CSV format:
% id, object_1, object_2, object_3, ...
%
% Each row is treated as a SET of objects.
% Repeated objects within the same row are counted only once for Jaccard.
% loop through categories

outputFilename = fullfile(pwd, '..', 'derivatives', 'group_level', 'jaccard_similarity.mat');

if exist(outputFilename, 'file') && cfg.force_recompute == false
    load(outputFilename)
else

    for category = cfg.categories
        category = char(category);

        for analysis_name = cfg.analysis_names
            analysis_name = char(analysis_name);

            %% Choose input file
            if strcmp(analysis_name, 'typical')
                filename = fullfile(pwd, '..', 'drawings', [category, '_objects.csv']);
            elseif strcmp(analysis_name, 'control')
                filename = fullfile(pwd, '..', 'drawings', [category, '_control_objects.csv']);
            elseif strcmp(analysis_name, 'photos')
                filename = fullfile(pwd, '..', 'photos', [category, '_objects.csv']);
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

            for i = 1:nRows
                rowObjects = objectData(i, :);
                rowObjects = strtrim(rowObjects);
                rowObjects = rowObjects(~ismissing(rowObjects) & rowObjects ~= "");

                % Treat as a set: remove duplicates within row
                rowObjects = unique(rowObjects);

                [tf, loc] = ismember(rowObjects, uniqueObjects);
                presence(i, loc(tf)) = true;
            end

            %% Compute pairwise Jaccard similarity
            jaccardSim_mat = zeros(nRows, nRows);

            for i = 1:nRows
                for j = 1:nRows
                    intersectionCount = sum(presence(i, :) & presence(j, :));
                    unionCount = sum(presence(i, :) | presence(j, :));

                    if unionCount == 0
                        jaccardSim_mat(i, j) = NaN;
                    else
                        jaccardSim_mat(i, j) = intersectionCount / unionCount;
                    end
                end
            end
            
            if cfg.dissimilarity
                jaccardSim_mat = 1 - jaccardSim_mat;
            end 

            jaccardSim.(analysis_name).(category).subject_mean.name = 'jaccardSim';
            jaccardSim.(analysis_name).(category).subject_mean.color = [0,0,0];
            jaccardSim.(analysis_name).(category).subject_mean.RDM = jaccardSim_mat;
            jaccardSim.(analysis_name).(category).subjects = rowIDs;

            disp("Computed pairwise Jaccard similarity matrix for:");
            disp([category, ' - ', analysis_name]);

        end
    end

    % Save results
    save(outputFilename, "jaccardSim")
end

% Write to results struct
d.DNN.jaccardSim = jaccardSim; % pretend this is a DNN approach

disp("Done.");

end