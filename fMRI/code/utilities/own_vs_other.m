function d = own_vs_other(d, cfg)

if ~isfield(cfg, 'fMRI_data_type');  cfg.fMRI_data_type = 'betas';end
cfg.correlation_type = 'Spearman';

% prepare
warning('off')
nmasks = numel(cfg.rois_of_interest);
nPred = length(cfg.analysis_names);
own_vs_other_mat = zeros(nPred, cfg.n, nmasks, length(cfg.categories));
own_cor = zeros(nPred, cfg.n, nmasks, length(cfg.categories), 1);
other_cor = zeros(nPred, cfg.n, nmasks, length(cfg.categories), cfg.n - 1);
var_names = {'brain_resp', 'stimulus', cfg.analysis_names{:}, 'own', 'ref_subject', 'pred_subject', 'category','roi'};
pred_formats = repmat({'double'}, 1, nPred);
big_fucking_table = table( ...
    'Size',[1, length(var_names)], ...
    'VariableTypes',{'double', 'string', pred_formats{:}, 'logical','string','string','string','string'}, ...
    'VariableNames',var_names);


figure;
tiledlayout(1, length(cfg.categories)+2)
title(['Own vs other - ', cfg.dnn])

%cycle through categories
for iCate = 1:length(cfg.categories)
    category = char(cfg.categories{iCate});

    % get all predictors
    allPredictors = zeros(nPred, cfg.nTrials/2, cfg.n);
    allPredictors_unique = allPredictors;
    for iPred = 1:nPred
        for iSub = 1:cfg.n
            allPredictors(iPred, :, iSub) = ...
                cell2mat({d.DNN.(cfg.dnn).stimuli.(category).(char(cfg.analysis_names{iPred}))(iSub).mean_sim});
        end
    end

    % residualized predictors
    if cfg.partial_cor && nPred > 1
        for iSub = 1:cfg.n

            sub_preds = squeeze(allPredictors(:, :, iSub))';
            uniquePred = zeros(size(sub_preds));

            for iPred = 1:nPred

                % all other predictors
                otherIdx = setdiff(1:nPred, iPred);
                X = sub_preds(:,otherIdx);
                X = [X, ones(size(X,1),1)];
                b = X \ sub_preds(:,iPred);

                % keep residual (unique variance)
                uniquePred(:,iPred) = sub_preds(:,iPred) - X*b;

            end
            allPredictors_unique(:, :, iSub) = uniquePred';
        end
    end

    all_mask_lables = cell(nmasks, 1);
    for j=1:nmasks

        % get mask name
        mask_label_short = cfg.rois_of_interest{j};
        all_mask_lables{j} = mask_label_short;

        % progress report
        disp(['Processing ROI: ', mask_label_short])

        for iSub = 1:length(cfg.subNums)

            % subject ID
            subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
            subID2 = sprintf('sub%0.3d', cfg.subNums(iSub));

            % get fMRI data for this subject
            stimNames = {d.betas.(category).(subID2).(mask_label_short).image_name}';
            if strcmp(cfg.fMRI_data_type, 'beta') % get mean beta response
                ownData = cell2mat({d.betas.(category).(subID2).(mask_label_short).mean_betas})';
            elseif strcmp(cfg.fMRI_data_type, 'pairDec') % get mean decoding data
                decRDM = d.pairRep.(subID2).(['w', mask_label_short]).rdm;
                decRDM(eye(size(decRDM))==1) = NaN;
                ownData = mean(decRDM, 2, 'omitnan');
                if iCate == 1
                    ownData = ownData(1:cfg.nTrials/2);
                elseif iCate == 2
                    ownData = ownData(cfg.nTrials/2+1:end);
                end
            end

            % add values to big tables
            current_row = height(big_fucking_table);
            rows = current_row+1:current_row+length(ownData);
            big_fucking_table.brain_resp(rows) = ownData;
            big_fucking_table.stimulus(rows) = stimNames;
            big_fucking_table.own(rows) = true;
            big_fucking_table.ref_subject(rows) = subID;
            big_fucking_table.pred_subject(rows) = subID;
            big_fucking_table.category(rows) = category;
            big_fucking_table.roi(rows) = mask_label_short;

            % correlate with stimuli features for unique predictor
            for iPred = 1:nPred

                % add to table
                big_fucking_table.(cfg.analysis_names{iPred})(rows) = squeeze(allPredictors(iPred, :, iSub))';

                % run correlation
                cor_vals = corr(ownData, squeeze(allPredictors_unique(iPred, :, :)), 'row', 'pairwise', 'type', cfg.correlation_type)';
                own_cor(iPred, iSub, j, iCate) = cor_vals(iSub);
                other_cor(iPred, iSub, j, iCate, :) = cor_vals(setdiff(1:cfg.n, iSub));
            end

            % compare own vs other
            own_vs_other_mat(:, iSub, j, iCate) = own_cor(:, iSub, j, iCate) - ...
                mean(squeeze(other_cor(:, iSub, j, iCate, :)), 2, 'omitnan');

            % write other subjects data to table
            for iSub_other = 1:length(cfg.subNums)
                if iSub == iSub_other
                    continue
                end

                current_row = height(big_fucking_table);
                rows = current_row+1:current_row+length(ownData);
                big_fucking_table.brain_resp(rows) = ownData;
                big_fucking_table.stimulus(rows) = stimNames;
                big_fucking_table.own(rows) = false;
                big_fucking_table.ref_subject(rows) = subID;
                big_fucking_table.pred_subject(rows) = sprintf('sub-%0.3d', cfg.subNums(iSub_other));
                big_fucking_table.category(rows) = category;
                big_fucking_table.roi(rows) = mask_label_short;

                for iPred = 1:nPred
                    big_fucking_table.(cfg.analysis_names{iPred})(rows) = squeeze(allPredictors(iPred, :, iSub_other))';
                end
            end
        end
    end

    % Mean and SEM
    X = squeeze(own_vs_other_mat(1, :, :, iCate));
    colMean = mean(X, 1, 'omitnan');
    rowsEM  = std(X, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(X),1));

    nexttile

    hold on;

    % Bar plot
    b = bar(colMean, 'FaceColor', [0.7 0.7 0.7]);

    % Error bars
    errorbar(1:size(X,2), colMean, rowsEM, ...
        'k', 'LineStyle', 'none', 'LineWidth', 1.5);

    % Individual data points
    for m = 1:size(X,2)

        % jitter x positions slightly
        xj = m + 0.15*(rand(size(X,1),1)-0.5);

        scatter(xj, X(:,m), ...
            40, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.6);

    end

    xlabel('ROI');
    xlim([0.5,nmasks + 0.5])
    ylabel('Own - other');
    box off;
    set(gca,'LineWidth',1.5,'FontSize',12);
    xticks(1:nmasks)
    xticklabels(all_mask_lables)
    title(category)
    hold off

end

% clean table
big_fucking_table = big_fucking_table(~any(ismissing(big_fucking_table), 2), :);

% LPFC_table = big_fucking_table(strcmp(big_fucking_table.roi, 'LPFC'), :);
% lme = fitlme(LPFC_table,...
%     ['brain_resp ~ typical + control + photos + own + typical:own + control:own + photos:own +' ...
%     '(typical + control + photos + own + typical:own + control:own + photos:own | ref_subject) + ' ...
%     '(typical + control + photos + own + typical:own + control:own + photos:own | pred_subject) + ' ...
%     '(typical + control + photos + own + typical:own + control:own + photos:own | category) + ' ...
%     '(typical + control + photos + own + typical:own + control:own + photos:own | stimulus)']);
% disp(lme)


% Mean and SEM
X = squeeze(mean(own_vs_other_mat(1, :, :, :), 4));
colMean = mean(X, 1, 'omitnan');
rowsEM  = std(X, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(X),1));

nexttile
hold on;

% Bar plot
b = bar(colMean, 'FaceColor', [0.7 0.7 0.7]);

% Error bars
errorbar(1:size(X,2), colMean, rowsEM, ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Individual data points
for m = 1:size(X,2)
    % jitter x positions slightly
    xj = m + 0.15*(rand(size(X,1),1)-0.5);
    scatter(xj, X(:,m), ...
        40, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.6);
end

%% Significance markers
yRange = range(ylim);
pvals = nan(1,nmasks);
tvals = nan(1,nmasks);

for m = 1:nmasks
    [~,p,~,stats] = ttest(X(:,m),0);
    pvals(m) = p;
    tvals(m) = stats.tstat;
end
[~,~,~,pvals_fdr] = fdr_bh(pvals);
for m = 1:nmasks
    if pvals_fdr(m) < 0.05
        % asterisk position
        yPos = colMean(m) + rowsEM(m) + 0.05*yRange;
        text(m, yPos, '*', ...
            'HorizontalAlignment','center', ...
            'FontWeight','bold', ...
            'FontSize',18);
    end
end

xlabel('ROI');
xlim([0.5,nmasks + 0.5])
ylabel('Own - other');
box off;
set(gca,'LineWidth',1.5,'FontSize',12);
xticks(1:nmasks)
xticklabels(all_mask_lables)
title('Catgeory average')
hold off

%% check correlation
own_cor_mean = mean(own_cor, 4);
X = squeeze(own_cor_mean(1, :, :));

% Mean and SEM
colMean = mean(X, 1, 'omitnan');
rowsEM  = std(X, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(X),1));

nexttile
hold on;

% Bar plot
b = bar(colMean, 'FaceColor', [0.7 0.7 0.7]);

% Error bars
errorbar(1:size(X,2), colMean, rowsEM, ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Individual data points
for m = 1:size(X,2)
    % jitter x positions slightly
    xj = m + 0.15*(rand(size(X,1),1)-0.5);
    scatter(xj, X(:,m), ...
        40, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.6);
end

xlabel('ROI');
xlim([0.5,nmasks + 0.5])
ylabel('Correlation');
box off;
set(gca,'LineWidth',1.5,'FontSize',12);
xticks(1:nmasks)
xticklabels(all_mask_lables)
title(['Correlation (mean ', cfg.fMRI_data_type ,'beta to sitmulus similarity)'])
hold off

warning('on')

end