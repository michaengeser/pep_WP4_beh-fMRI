function permRes = doPermutations(RDMs, random_seqs, cfg)

% define configurations
rng('default')
rng(1)
cfg.plotting = false;
cfg.plot_rdm = false;
permRes = zeros(numel(RDMs), numel(RDMs), cfg.n_permutations);

for iRDM = 1:numel(RDMs)

    permutation_RDMs = RDMs(iRDM:end);

    for perm = 1:cfg.n_permutations

        % randomize RDM
        if strcmp(cfg.permutation_type, 'row_col_shuffle_ref')
            % shuffle the order of rows and columns
            ref_RDM = permutation_RDMs(1);
            ref_RDM.RDM = ref_RDM.RDM(random_seqs{perm}, random_seqs{perm});
            % replace reference RDM by permutated RDM
            permutation_RDMs(1) = ref_RDM;
        else
            error('Permutation type curretly not supported, use row_col_shuffle_ref')
        end


        % get correaltion of RDMs matrix
        if cfg.partial_cor
            % partial correlaltion
            [~, res_rdm, ~ , cfg] = partial_cor_RDM(cfg, permutation_RDMs);
        else
            % regular correlation
            [~, res_rdm, ~] = cor_RDM(permutation_RDMs, cfg);
        end

        % store in output rdm
        permRes(iRDM:end, iRDM, perm) = res_rdm(:, 1);

        if mod(perm/cfg.n_permutations, 0.1) == 0
            disp([num2str((perm/cfg.n_permutations)*100), '% of permutations of ',...
                RDMs(iRDM).name, ' are done'])
        end

    end
end
end