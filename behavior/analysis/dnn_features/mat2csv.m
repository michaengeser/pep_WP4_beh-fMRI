% Script to convert .mat files (struct with strings + numbers) into CSVs

% Root folder to search
rootFolder = pwd;

% Recursively find all .mat files
matFiles = dir(fullfile(rootFolder, '**', 'exp_1', '*.mat'));

for k = 1:length(matFiles)
    % Full path to .mat file
    matFileName = matFiles(k).name;
    matPath = fullfile(matFiles(k).folder, matFileName);

    % Load the .mat file (assumes it contains a single struct)
    data = load(matPath);
    data = data.dnn_features.all_net_act;

    % Extract struct (assume only one variable inside)
    fn = fieldnames(data);
    s = data.(fn{1});

    T = table;
    for i = 1:width(data)
        T.(data(i).image_name) = data(i).net_act;
    end


    % Generate output filename using keywords from path
    keyword = [];
    relPath = erase(matPath, rootFolder);
    relParts = strsplit(relPath, filesep);
    relParts(cellfun('isempty', relParts)) = []; % remove empties
    keyword{1} = relParts{1};

    if contains(matFileName, 'Original', 'IgnoreCase', true)
        keyword{2} = 'original';
    else
        keyword{2} = 'draw3D';
    end

    if contains(matFileName, 'own', 'IgnoreCase', true)
        keyword{3} = 'own';
    elseif contains(matFileName, 'control', 'IgnoreCase', true)
        keyword{3} = 'control';
    elseif contains(matFileName, 'photo', 'IgnoreCase', true)
        keyword{3} = 'photo';
    end


    if contains(matFileName, 'bathroom', 'IgnoreCase', true)
        keyword{4} = 'bathroom';
    elseif contains(matFileName, 'kitchen', 'IgnoreCase', true)
        keyword{4} = 'kitchen';
    end

    % Compose file name
    outName = strjoin(keyword, '_');

    % Save as CSV in same folder (or change destination if needed)
    outPath = fullfile(pwd, [outName '.csv']);
    writetable(T, outPath);

    fprintf('Saved CSV: %s\n', outPath);
end
