function [d,cfg] = stimuli_predictor_object_overlap(d,cfg)

%% Defaults
if ~isfield(cfg,'categories')
    cfg.categories = {'kitchen','bathroom'};
end

if ~isfield(cfg,'analysis_names')
    cfg.analysis_names = {'typical','control'};
end

%% Loop categories

for iCate = 1:numel(cfg.categories)

    category = cfg.categories{iCate};

    % Load stimulus object file
    stimFile = fullfile(pwd, '..', 'stimuli', category,...
        [category '_objects.csv']);

    opts.Delimiter = ",";
    opts = delimitedTextImportOptions("NumVariables", 100);
    stimTab = readtable(stimFile, opts);

    stimNames = stimTab{2:end,1};
    stimObjects = stimTab{2:end,2:end};

    % Loop predictor types
    for iAnalysis = 1:numel(cfg.analysis_names)

        analysis_name = cfg.analysis_names{iAnalysis};

        % Predictor file
        if strcmp(analysis_name,'typical')
            predFile = fullfile(pwd, '..', 'drawings', [category, '_objects.csv']);
        elseif strcmp(analysis_name, 'control')
            predFile = fullfile(pwd, '..', 'drawings', [category, '_control_objects.csv']);
        elseif strcmp(analysis_name, 'photos')
            predFile = fullfile(pwd, '..', 'photos', [category, '_objects.csv']);
        else
            error('Unknown analysis type')
        end

        predTab = readtable(predFile,'TextType','string', 'Format','auto');

        predNames = predTab{:,1};
        predObjects = predTab{:,2:end};

        % Subject loop

        for sub = 1:cfg.n

            subNum = cfg.subNums(sub);

            %% find subject predictor image

            predIdx = predNames == subNum;

            if ~any(predIdx)
                warning('No predictor image found for sub-%03d',subNum)
                continue
            end

            currentPredObjects = predObjects(predIdx, :);
            currentPredObjects = currentPredObjects(~strcmp((currentPredObjects), ""));

            % Compare with all stimuli
            overlapVals = nan(height(stimObjects),1);

            for iStim = 1:height(stimObjects)

                % get stimulus object
                stimSet = stimObjects(iStim, :);
                stimSet = stimSet(~ismissing(stimSet));
                stimSet = stimSet(~strcmp((stimSet), ""));

                % get overlap
                interSize = numel(intersect(stimSet,currentPredObjects));
                unionSize = numel(union(stimSet,currentPredObjects));

                if unionSize==0
                    overlapVals(iStim)=NaN;
                else
                    overlapVals(iStim)=interSize/unionSize;
                end

            end

            %% store
            d.DNN.objectSim.stimuli.(category).(char(analysis_name))(sub).mean_sim = ...
                overlapVals;
            d.DNN.objectSim.stimuli.(category).(char(analysis_name))(sub).image_names = ...
                stimNames;
        end

        disp(['Finished ' analysis_name ' - ' category])

    end
end

end

