function d = compare_task_RDMs_to_predictor_RDMs(d, cfg)
%  COMPARE_TASK_RDMS_TO_PREDICTOR_RDMS Brief summary of this function.
%
% Detailed explanation of this function.
% evaluate input
if ~isfield(cfg, 'predictor_RDMs'); cfg.predictor_RDMs = {'Typcial Drawings Style', 'Photos vgg16_imagenet Late',...
        'Control Images vgg16_imagenet Late', 'Typical Images vgg16_imagenet Late', 'Survey responses'}; end
if ~isfield(cfg, 'partial_cor'); cfg.partial_cor = true;end
if cfg.partial_cor
    if ~isfield(cfg, 'RDM_to_partial_out'); cfg.RDM_to_partial_out = cfg.predictor_RDMs; end
end
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'Pearson';end
if ~isfield(cfg, 'plot_rdm'); cfg.plot_rdm = false;end
if ~isfield(cfg, 'permutation_test'); cfg.permutation_test = false;end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 10000;end
if ~isfield(cfg, 'bootstrapping'); cfg.bootstrapping = false;end
if ~isfield(cfg, 'n_bootstrapp_iterations'); cfg.n_bootstrapp_iterations = 10000;end
if ~isfield(cfg, 'add_legend'); cfg.add_legend = false;end
if ~isfield(cfg, 'show_single_cate'); cfg.show_single_cate = false;end
if ~isfield(cfg, 'permutation_type'); cfg.permutation_type = 'row_col_shuffle_ref';end
if ~isfield(cfg, 'bootstrapp_type'); cfg.bootstrapp_type = 'removing';end
if ~isfield(cfg, 'order_predictors'); cfg.order_predictors = false;end
if ~isfield(cfg, 'save_name'); cfg.save_name = 'compare_task_RDMs_to_predictor_RDMs';end
if ~isfield(cfg, 'xaxis_labels'); cfg.xaxis_labels = true;end
if ~isfield(cfg, 'plot_type'); cfg.plot_type = 'violin';end
if ~isfield(cfg, 'scatter_in_violin'); cfg.scatter_in_violin = 0;end
if cfg.scatter_in_violin == 0
    violin_type = 'full';
elseif cfg.scatter_in_violin == 1
    violin_type = 'half';
end

if ~isfield(cfg, 'task_plotting'); cfg.task_plotting = true;end
if ~isfield(cfg, 'task_plotting'); cfg.plott_gap = 0;end
if ~isfield(cfg, 'plotting_predictors'); cfg.plotting_predictors = 1:numel(cfg.predictor_RDMs);end
cfg.plot_rdm = false;
% if permutations with sign flips are used, corrleation type needs to be
% spearman
if strcmp(cfg.permutation_type, 'sign_flip_ref')
    cfg.partial_correlation_type = 'spearman';
end
% get variable attributes
colors = zeros(numel(cfg.RDM_to_partial_out),3);
short_names = cell(1, numel(cfg.RDM_to_partial_out));
for var = 1:numel(cfg.RDM_to_partial_out)
    if strcmp(cfg.RDM_to_partial_out{var}, 'Typical Images vgg16_imagenet Late')
        short_names{var} = 'Typcial drawing';
        colors(var,:) = [1, 0, 1];
    elseif strcmp(cfg.RDM_to_partial_out{var}, 'Control Images vgg16_imagenet Late')
        short_names{var} = 'Control drawing';
        colors(var,:) = [.8, .8, .8];
    elseif strcmp(cfg.RDM_to_partial_out{var}, 'Typcial Drawings Style')
        short_names{var} = 'Drawing style ratings';
        colors(var,:) = [.5, .5, .5];
    elseif strcmp(cfg.RDM_to_partial_out{var}, 'Typcial Drawings Content')
        short_names{var} = 'Typical drawing ratings';
        colors(var,:) = [.7, 0, .7];
    elseif strcmp(cfg.RDM_to_partial_out{var}, 'Photos vgg16_imagenet Late')
        short_names{var} = 'Photos';
        colors(var,:) = [.4, .9, 1];
    end
end
% prepare tiled figure
figure;
hold on
previous_x_pos = 0;
% tiledlayout(numel(cfg.tasks_of_interest)/2,2)
% prepare random permutation and/or bootstrapping (each task and category should have the same
% random samplings)
if cfg.permutation_test
    % generate permutated subjects list
    rng(0)
    if cfg.partial_cor
        random_seqs = cell(numel(cfg.RDM_to_partial_out), cfg.n_permutations);
    else
        random_seqs = cell(1, cfg.n_permutations);
    end
    for i = 1:cfg.n_permutations
        if ismember(cfg.permutation_type, {'row_col_shuffle_ref',...
                'row_col_shuffle_pred', 'row_col_shuffle_pred_all', 'row_col_shuffle_pred_plus_reoder'})
            % get random sequence of rows and columns
            if cfg.partial_cor
                for ii = 1:numel(cfg.RDM_to_partial_out)
                    random_seqs{ii,i} = randperm(cfg.n);
                end
            else
                random_seqs{i} = randperm(cfg.n);
            end
        elseif strcmp(cfg.permutation_type, 'sign_flip_ref')
            % get random rows and columns to flip sign
            random_seqs{1,i} = randperm(cfg.n, randi([1, cfg.n]));
        end
    end
end

% loop through tasks
for voi_n = 1:numel(cfg.tasks_of_interest)
    voi = char(cfg.tasks_of_interest(voi_n));
    % loop through categories
    for cate_num = 1:numel(cfg.categories)
        category = char(cfg.categories{cate_num});
        % get task RDM names
        all_ref_names = {d.([category,'_RDM']).ratingRDM.name};
        % get voi RDM
        ref_idx = find(strcmp(all_ref_names, voi));
        if ~isempty(ref_idx)
            ref_RDM = d.([category,'_RDM']).ratingRDM(ref_idx);
        else
            tem_preditor_RDMs = cfg.predictor_RDMs; % store predictors temporally
            cfg.predictor_RDMs = {voi};
            RDMs = d.([category,'_RDM']).ratingRDM(1); % fill with some RDM as placeholder
            labels = {RDMs.name};
            [ref_RDM, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, category);
            cfg.predictor_RDMs = tem_preditor_RDMs; % write back predictors
            ref_RDM = ref_RDM(2:end); % remove placeholder
            cfg.labels = cfg.labels{2:end};
        end
        % get canditate/predictor RDMs
        RDMs = ref_RDM;
        labels = {ref_RDM.name};
        [RDMs, cfg.labels] =  evaluate_predictor_RDMs(d, RDMs, labels, cfg, category);
        % make a cell that holds all predictor RDM structs
        for field = 1:numel({RDMs.name})
            RDMs(field).name = char(cfg.labels{field}); % give it a comprehensive name
        end
        % partial correlation
        if cfg.partial_cor
            [~, r_mat, ~, cfg] = partial_cor_RDM(cfg, RDMs);
        else
            [~, r_mat, ~] = cor_RDM(RDMs,cfg);
        end
        % store results in table
        res_table = table;
        res_table.name = cfg.labels(2:end)';
        res_table.r_val = r_mat(2:end, 1);
        % store in data struct
        d.compare_task_to_predictor.(voi).(category) = res_table;

        if cfg.permutation_test
            nPred = numel(cfg.RDM_to_partial_out);

            % reference RDM
            ref_RDM = RDMs(1).RDM;

            % get predictors and regress out from each other
            for iPred = 1:nPred
                RDMs(iPred+1).RDM(eye(cfg.n)==1) = 0;
                all_preds(:,iPred) = squareform(RDMs(iPred+1).RDM);
            end

            % residualized predictors
            uniquePred = zeros(size(all_preds));
            if cfg.partial_cor && nPred > 1
                for iPred = 1:nPred

                    % all other predictors
                    otherIdx = setdiff(1:nPred, iPred);
                    X = all_preds(:,otherIdx);
                    X = [X, ones(size(X,1),1)];
                    b = X \ all_preds(:,iPred);

                    % keep residual (unique variance)
                    uniquePred(:,iPred) = all_preds(:,iPred) - X*b;

                end
            else
                uniquePred = all_preds;
            end

            % get actual correaltion
            obs_r = r_mat(2:end, 1);

            fprintf('Running %d permutations...\n',cfg.n_permutations)

            % create permutations
            allPerms = nan(cfg.n_permutations, nchoosek(cfg.n, 2));
            parfor p = 1:cfg.n_permutations

                % shuffle neural RDM
                RDM_shuffled = ref_RDM(random_seqs{1,p}, random_seqs{1,p});
                allPerms(p, :) = squareform(RDM_shuffled);
            end

            % run correlations for each predictor
            d.compare_task_to_predictor.permutation_test.(voi).(category) = ...
                corr(allPerms', uniquePred, 'row', 'pairwise', 'type', cfg.partial_correlation_type)';
        end

        % runtime control
        disp(['Compare inter-subject RDM of ', voi, ' with correaltion of predictor RDMs - ', category])

        % store in data struct
        d.compare_task_to_predictor.(voi).(category) = res_table;
    end

    % average categories
    res_table.r_val_cate1 = d.compare_task_to_predictor.(voi).(cfg.categories{1}).r_val;
    res_table.r_val_cate2 = d.compare_task_to_predictor.(voi).(cfg.categories{2}).r_val;
    res_table.r_val = (d.compare_task_to_predictor.(voi).(cfg.categories{1}).r_val +...
        d.compare_task_to_predictor.(voi).(cfg.categories{2}).r_val)/2;
    d.compare_task_to_predictor.(voi).category_average = res_table;

    % get confidence intervals
    if cfg.permutation_test
        % get p values of random permutation
        perm_r_mat = (d.compare_task_to_predictor.permutation_test.(voi).(cfg.categories{1}) +...
            d.compare_task_to_predictor.permutation_test.(voi).(cfg.categories{2}))/2;
        d.compare_task_to_predictor.permutation_test.(voi).category_average = perm_r_mat;
    end


    % plotting

    % filter predictors to plot
    if numel(cfg.plotting_predictors) ~= numel(cfg.predictor_RDMs)
        res_table = res_table(cfg.plotting_predictors, :);
    end

    if ~isempty(short_names) && numel(cfg.plotting_predictors) ~= numel(short_names)
        short_names = short_names(cfg.plotting_predictors);
        colors = colors(cfg.plotting_predictors, :);
    end

    if numel(cfg.plotting_predictors) ~= height(perm_r_mat)
        perm_r_mat = perm_r_mat(cfg.plotting_predictors, :);
    end

    % loop through res_table and add according variables
    if cfg.plotting_predictors == 1
        if strcmp(voi, 'typicality')
            res_table.color_R = .98;
            res_table.color_G = .41;
            res_table.color_B = .91;
        elseif strcmp(voi, 'familiarity')
            res_table.color_R = .87;
            res_table.color_G = .01;
            res_table.color_B = .87;
        elseif strcmp(voi, 'aesthetics')
            res_table.color_R = .6;
            res_table.color_G = .02;
            res_table.color_B = .6;
        elseif strcmp(voi, 'usability')
            res_table.color_R = .47;
            res_table.color_G = .03;
            res_table.color_B = .45;
        elseif strcmp(voi, 'complexity')
            res_table.color_R = .24;
            res_table.color_G = .04;
            res_table.color_B = .3;
        elseif strcmp(voi, 'IES')
            res_table.color_R = 1;
            res_table.color_G = 0;
            res_table.color_B = 1;
        end

    else
        clr = summer(numel(RDMs));
        for row = 1:height(res_table)
            lowerCaseName = lower(res_table.name{row});
            % colors
            if ~isempty(colors)
                res_table.color_R(row) = colors(row,1);
                res_table.color_G(row) = colors(row,2);
                res_table.color_B(row) = colors(row,3);
            else % if no colors specified make the bars in colomap defined in clr
                if contains(lowerCaseName, {'control', 'typical', 'own'})
                    if contains(lowerCaseName, {'typical', 'own'})
                        res_table.color_R(row) = 1;
                        res_table.color_G(row) = 0;
                        res_table.color_B(row) = 1;
                    elseif contains(lowerCaseName, {'control'})
                        res_table.color_R(row) = .8;
                        res_table.color_G(row) = .8;
                        res_table.color_B(row) = .8;
                    end
                else
                    res_table.color_R(row) = clr(row,1);
                    res_table.color_G(row) = clr(row,2);
                    res_table.color_B(row) = clr(row,3);
                end
            end
        end
    end
    % short names
    for row = 1:height(res_table)
        if ~isempty(short_names)
            res_table.short_names{row} = short_names{row};
        else
            if contains(lowerCaseName, cfg.dnns)
                for i = 1:length(cfg.dnns)
                    if contains(lowerCaseName, cfg.dnns{i})
                        res_table.short_names{row} = strrep(cfg.dnns{i}, '_', ' ');
                    end
                end
            else
                res_table.short_names{row} = cfg.labels{row+1};
            end
        end

        % p val
        if cfg.permutation_test
            % get p value from random permutation (one-sided test aginst
            % permutation distribution)
            p_value = sum(perm_r_mat(row,:) >= res_table.r_val(row)) / cfg.n_permutations;
            res_table.p_val(row) = p_value;
            % get p values for categories
            res_table.p_val_cate1(row) = sum(d.compare_task_to_predictor.permutation_test.(voi).(cfg.categories{1})(row,:)...
                >= res_table.r_val_cate1(row)) / cfg.n_permutations;
            res_table.p_val_cate2(row) = sum(d.compare_task_to_predictor.permutation_test.(voi).(cfg.categories{2})(row,:)...
                >= res_table.r_val_cate2(row)) / cfg.n_permutations;
            res_table.ci_upper(row) = res_table.r_val(row) - prctile(perm_r_mat(row,:), 5);
            res_table.ci_lower(row) = res_table.r_val(row) - prctile(perm_r_mat(row,:), 95);
        elseif cfg.bootstrapping
            % get p value from randomly sampled data (one-sided test of
            % bootstrapping distribution against 0)
            p_value = sum(bootstrapping_r_vals(row,:) <= 0) / cfg.n_bootstrapp_iterations;
            res_table.p_val(row) = p_value;
            % get p values for categories
            res_table.p_val_cate1(row) = sum(boot_r_vals_cate1(row,:) <= 0) / cfg.n_bootstrapp_iterations;
            res_table.p_val_cate2(row) = sum(boot_r_vals_cate2(row,:) <= 0) / cfg.n_bootstrapp_iterations;
        else
            % get p values from r values
            N = nchoosek(cfg.n, 2);
            res_table.p_val(row) = r2p(res_table.r_val(row), N);
            res_table.p_val_cate1(row) = r2p(res_table.r_val_cate1(row), N);
            res_table.p_val_cate2(row) = r2p(res_table.r_val_cate2(row), N);
        end
    end


    % order row based on r values
    if cfg.order_predictors
        res_table = sortrows(res_table, 'r_val', 'descend');
    end

    % store in data struct
    d.compare_task_to_predictor.(voi).category_average = res_table;

    % make plot
    if strcmp(cfg.plot_type, 'violin')

        for xiPos = 1:height(res_table)
            current_x_pos = previous_x_pos + xiPos;
            if cfg.permutation_test
                Y = {perm_r_mat(xiPos,:)'};
            end
            % make violin plot

            % make violin plot
            currentColor = [res_table.color_R(xiPos), res_table.color_G(xiPos), res_table.color_B(xiPos)];
            mainHandles(current_x_pos) = daviolinplot(Y,...
                'color', currentColor,...
                'violin', violin_type, 'violinalpha', 0.2,...
                'scatter',cfg.scatter_in_violin,'scatteralpha',0.2,'jitter',1,'scattercolors', 'same', 'scattersize', 5,...
                'box', 0,...
                'outliers',0);
            mainHandles(current_x_pos).ds.Vertices(:, 1) = mainHandles(current_x_pos).ds.Vertices(:, 1) + current_x_pos - 1;
            mainHandles(current_x_pos).ds.EdgeColor = currentColor;
            mainHandles(current_x_pos).ds.EdgeAlpha = 0.5;
            mainHandles(current_x_pos).ds.LineWidth = 2;
            if cfg.scatter_in_violin == 1
                mainHandles(current_x_pos).sc.XData = mainHandles(current_x_pos).sc.XData + current_x_pos - 1.05;
                mainHandles(current_x_pos).sc.MarkerEdgeColor = currentColor;
                mainHandles(current_x_pos).sc.MarkerEdgeAlpha = 0.2;
            end
        end


    elseif strcmp(cfg.plot_type, 'bar')
        for xiPos = 1:height(res_table)
            current_x_pos = previous_x_pos + xiPos;
            % Draw individual bar
            barColor = [res_table.color_R(xiPos),res_table.color_G(xiPos),res_table.color_B(xiPos)];
            mainHandles(current_x_pos) = bar(current_x_pos, res_table.r_val(xiPos), 'FaceColor', barColor, 'EdgeColor', 'k');
        end
    end
    for xiPos = 1:height(res_table)
        current_x_pos = previous_x_pos + xiPos;
        if strcmp(cfg.plot_type, 'bar')
            if contains(string(res_table.name(xiPos)), 'Original')
                % Apply hatch only to this bar
                hatchfill2(mainHandles(current_x_pos), 'HatchAngle', 45, ...
                    'HatchColor', 'k', ...
                    'HatchLineWidth', 1);
            end
        end
        if cfg.permutation_test
            if strcmp(cfg.plot_type, 'bar')
                % add confidence interval if available
                if ismember('ci_lower', res_table.Properties.VariableNames)
                    r_val = res_table.r_val(xiPos);
                    errorHandles(current_x_pos) = errorbar(current_x_pos, r_val,...
                        r_val-res_table.ci_lower(xiPos), res_table.ci_upper(xiPos)-r_val, 'k', 'LineWidth', 1);  % Error bars
                end
                % add bootstrapping median if available
                if ismember('boot_median', res_table.Properties.VariableNames)
                    text(current_x_pos, res_table.boot_median(xiPos), '-', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
                end
            elseif strcmp(cfg.plot_type, 'violin')
                % plot oberseved mean r
                plot([current_x_pos-0.15,current_x_pos+0.15], [res_table.r_val(xiPos), res_table.r_val(xiPos)],...
                    'Color', [res_table.color_R(xiPos), res_table.color_G(xiPos), res_table.color_B(xiPos)], 'LineWidth',3);
            end
        end
        % add marks for single category
        if cfg.show_single_cate
            if cfg.exp_num == 1
                cate_mark1 = 'B';
                cate_mark2 = 'K';
            elseif cfg.exp_num == 2
                cate_mark1 = 'B';
                cate_mark2 = 'L';
            end
            text(current_x_pos-0.2, res_table.r_val_cate1(xiPos), cate_mark1, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 5);
            text(current_x_pos-0.2, res_table.r_val_cate2(xiPos), cate_mark2, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 5);
        end

    end
    % make gap between reference RDMs
    previous_x_pos = current_x_pos + cfg.plott_gap;
end

if cfg.plotting

    % collect p values
    all_p_vals = nan(height(res_table), numel(cfg.tasks_of_interest));
    for voi_i = 1:numel(cfg.tasks_of_interest)
        voi = char(cfg.tasks_of_interest(voi_i));
        all_p_vals(:, voi_i) = d.compare_task_to_predictor.(voi).category_average.p_val;
    end

    % get asterisks
    all_asterisks = cell(height(res_table), numel(cfg.tasks_of_interest));
    for i_pred = 1:height(res_table)
        % do fdr correction
        [~, ~, ~, fdr_pval] = fdr_bh(all_p_vals(i_pred, :));
        all_asterisks(i_pred, :) = pval2asterisks(fdr_pval, 'none');

        % write back adjusted p values and print them
        disp([newline, newline])
        disp(string(res_table.name(i_pred)))
        for voi_i = 1:numel(cfg.tasks_of_interest)
            voi = char(cfg.tasks_of_interest(voi_i));
            d.compare_task_to_predictor.(voi).category_average.p_val(i_pred) = fdr_pval(voi_i);
            disp(['FDR corrected p value for ', voi, ': ', num2str(fdr_pval(voi_i)),...
                ' ', char(all_asterisks(i_pred, voi_i))])
            disp(['R value for ', voi, ': ', ...
                num2str(d.compare_task_to_predictor.(voi).category_average.r_val((xiPos)))])
        end
    end

    % plot asterisks
    ast_vec = reshape(all_asterisks, 1, []);
    y_pos = zeros(1, length(mainHandles));
    count = 0;
    for i_bar = 1:length(mainHandles)
        if strcmp(cfg.plot_type, 'bar')
            if ~isgraphics(mainHandles(i_bar))
                continue
            end

            % get y position based on error bars
            y_pos(i_bar) = errorHandles(i_bar).YPositiveDelta + errorHandles(i_bar).YData + 0.08;

        elseif strcmp(cfg.plot_type, 'violin')
            if ~isgraphics(mainHandles(i_bar).ds)
                continue
            end
            y_pos(i_bar) = max(mainHandles(i_bar).ds.Vertices(:, 2)) + 0.08;

        end
        count = count + 1;
        astHandle{count} = text(i_bar, y_pos(i_bar), ast_vec{count}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 20);
    end

    for i_ast = 1:length(mainHandles)
        astHandle{i_ast}.Position(2) = max(y_pos);
    end

    % get aesthetics
    hold off
    if cfg.partial_cor
        ylabel(['Partial correlation [r]', newline]);
    else
        ylabel([cfg.correlation_type, ' correlation [r]', newline]);
    end
    %title('Compare reference RDM with predictors')
    if isfield(cfg, 'plot_type')
        ylim(cfg.ylim)
    else
        ylim([-0.1, max(res_table.r_val) + 0.1])
    end
    xlim([0, previous_x_pos + 1])
    set(gca, 'LineWidth', 2, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold')
    ax = gca;
    ax.Box = 'off';
    yline(0, 'LineWidth', 2, 'Color', 'k');
    % get labels
    if cfg.xaxis_labels
        if cfg.task_plotting
            xticks(ceil(length(cfg.plotting_predictors)/2):...
                length(cfg.plotting_predictors)+cfg.plott_gap:...
                (length(cfg.plotting_predictors)+2)*length(cfg.tasks_of_interest));
            xticklabels(cfg.tasks_of_interest);
            xtickangle(45);
        else
            xticks(1:height(res_table));
            xticklabels(res_table.short_names);
            xtickangle(45);
        end
    else
        xticklabels([]);
        ax.XColor = 'none';
    end
    % add legend to last plot
    if cfg.add_legend
        legend(res_table.short_names, 'Location','northeastoutside');
    end
    % saving
    % fig_path = fullfile(pwd, 'figures', ['exp_', num2str(cfg.exp_num)], 'compare_task_RDMs_to_predictor_RDMs');
    % save_plot(cfg.save_name, fig_path)
end
end