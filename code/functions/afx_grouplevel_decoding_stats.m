% function [comp_results, mask_total] = afx_grouplevel_decoding_stats(mask_ACCs, subs_todo, m_ind)
%
% Computes t-test statistics for current ACC matrix.
%
% IN: All accuracy matrices for one mask (all_mask_ACCs), subjects to do
% (subs_todo), mask index nr (m_ind), targetfolder for saving. If present,
% will use chancelevel input to determine testing against chance.
% OUT: Result struct with all needed statistics for plotting.
%
% Needs: afx_timeresolved_decoding_grouplevel input
%
% Kai and Ingmar, Jun 2020

% Expanded by Ingmar, September 2020

function [comp_results, mask_total] = afx_grouplevel_decoding_stats(run_cfg, all_mask_ACCs, m_ind)
%% Info

n_tests         = length(all_mask_ACCs.ACCs{1}.axislabel.testtime_xaxis_vals); % # of tests for MCC
comp_results    = cell(size(all_mask_ACCs.ACCs)); % Different # of comparisons per mask
subs_todo       = run_cfg.subs_todo; % subs to process

%% Loop
% Loop here, within function, to better access ACC fields per comparison
% and add mask info to all comparisons

for a_ind = 1:length(all_mask_ACCs.ACCs) % class comparisons
    
    result = struct(); % init
    
    % Some stat settings
    if isfield(run_cfg.results, 'chancelevel') && ~isempty(run_cfg.results.chancelevel) && ...
            0 < run_cfg.results.chancelevel && run_cfg.results.chancelevel < 100
        chancelvl = run_cfg.results.chancelevel; % set
    elseif contains(all_mask_ACCs.ACCs{1}.title, 'vs')
        chancelvl = 100/(length(strfind(all_mask_ACCs.ACCs{1}.title, 'vs')) + 1); % set
    elseif contains(all_mask_ACCs.ACCs{1}.title, 'all', 'IgnoreCase', true)
        chancelvl = 100/3; % assume it's three-class
        warningv('afx_grouplevel_decoding_stats:threeclass_assumption', ...
            'afx_grouplevel_decoding_stats has no info on class comparisons, is assuming three-class decoding!');
    end
    base_alpha  = 0.05/2; % two-sided
    dispv(1, '    Base alpha is %f, corrected for 2-sided testing (above & below chance)', base_alpha);
    
    % ROI MCC correction (bonf)
    if strcmpi(run_cfg.data_type, 'ROI') && (~isfield(run_cfg.results, 'ROI_MCC') || ...
            isempty(run_cfg.results.ROI_MCC) || run_cfg.results.ROI_MCC == 1)
        base_alpha = base_alpha/numel(run_cfg.ROI_masks); % bonferroni
        dispv(1, ['    Bonferroni-correction for the number of ROI tests was ', ...
            'applied to the alpha value. To disable this, set run_cfg.results.ROI_MCC == 0'])
    elseif strcmpi(run_cfg.data_type, 'ROI')
        dispv(1, ['    Bonferroni-correction for the number of ROI tests was ', ...
            'not applied to the alpha because run_cfg.results.ROI_MCC == 0'])
    end
    
    % Names & titles
    curr_maskname = all_mask_ACCs.mask_name;
    group_resname = [curr_maskname '_comp-' num2str(a_ind) '_' all_mask_ACCs.ACCs{a_ind}.title]; % comparison
    group_resname = [group_resname, '_n-' num2str(length(subs_todo))]; % number of subjects
    result.group_resname = strrep(strrep(group_resname, ' ', '_'), '.', '_');
    fprintf('Getting stats for result.group_resname: %s (comp %i/%i)\n', result.group_resname, a_ind, length(all_mask_ACCs.ACCs))
    
    % Get current data, subjects as first dim
    result.ACC_elements = shiftdim(all_mask_ACCs.ACCs{a_ind}.ACC, 2);
    assert(~any(isnan(result.ACC_elements), 'all'), ...
        'There are NaNs present in result.ACC_elements, check why!');
    
    % Info
    result.createdate       = datestr(now);
    result.full_datafile    = all_mask_ACCs.full_datafile;
    result.mask_name        = curr_maskname;
    result.mask_index       = m_ind;
    result.comparison_ind   = a_ind;
    result.comparison_str   = sprintf('Comparison %i: %s', a_ind, all_mask_ACCs.ACCs{a_ind}.title);
    result.subs_todo        = subs_todo;
    result.grouplvlversion  = 'v Ingmar, 14-Aug-21';
    
    % Labels
    result.axislabels   = all_mask_ACCs.ACCs{a_ind}.axislabel;
    result.title        = all_mask_ACCs.ACCs{a_ind}.title;
    
    % Get measure of data distribution across subjects
    result.kurtosis_all    = mean(kurtosis(result.ACC_elements, 1, 1), 'all');
    result.skewness_all    = mean(skewness(result.ACC_elements, 1, 1), 'all');
    
    % Descriptives
    result.mean             = squeeze(mean(result.ACC_elements, 1));
    result.diag.mean        = diag(result.mean);
    [result.min, min_ind]   = min(result.mean, [], 'all', 'linear'); % Beware: min_ind/max_ind are linear indices
    [result.max, max_ind]   = max(result.mean, [], 'all', 'linear');
    result.chancelevel      = chancelvl;
    
    %% Uncorrected t-test stats
    
    dispv(1, '  A. Performing uncorrected %s-sided t-test for comparison...', 'right');
    [H_uncorr,P_uncorr,~, STATS] = ttest(result.ACC_elements, ...
        chancelvl, 'alpha', base_alpha, 'tail', 'right');
    [~,~,CItwosided,~]  = ttest(result.ACC_elements, chancelvl, 'alpha', ...
        base_alpha); % two-sided for CI
    
    % add sig. below chance accuracies as -1 to H, unless p-val > 0.95
    % is due to being extremely close to chance
    H_uncorr(P_uncorr>(1-base_alpha)) = -1;
    isnull = squeeze(round(P_uncorr, 4)) == 1 & squeeze(round(result.mean, 4)) == ...
        round(chancelvl, 4);
    H_uncorr(isnull == 1) = 0;
    
    % Treat NaNs (likely due to extremely similar or exactly the same
    % values) as zero
    if any(isnan(H_uncorr), 'all')
        H_uncorr(isnan(H_uncorr)) = 0;
    end
    
    result.uncorr.H            = squeeze(H_uncorr);
    result.uncorr.P            = squeeze(P_uncorr);
    result.uncorr.CItwosided   = squeeze(CItwosided);
    result.uncorr.CItwos_diff  = squeeze(diff(CItwosided));
    result.uncorr.alpha        = base_alpha;
    result.uncorr.SD           = squeeze(STATS.sd);
    result.uncorr.stderr       = result.uncorr.SD/sqrt(size(result.mean, 2));
    result.uncorr.DF           = squeeze(STATS.df(1));
    result.uncorr.T            = squeeze(STATS.tstat);
    result.uncorr.Tmin         = result.uncorr.T(min_ind);
    result.uncorr.Tmax         = result.uncorr.T(max_ind);
    result.uncorr.stderr_minval= result.uncorr.stderr(min_ind);
    result.uncorr.stderr_maxval= result.uncorr.stderr(max_ind);
    
    % Stats relating to the min/max means
    % Beware: min_ind/max_ind are linear indices
    result.uncorr.CItwosided_min_meanval   = CItwosided(:, min_ind);
    result.uncorr.CItwosided_max_meanval   = CItwosided(:, max_ind);
    result.uncorr.CItwosided_max_diff      = max(result.uncorr.CItwos_diff);
    result.uncorr.CItwosided_min_diff      = min(result.uncorr.CItwos_diff);
    
    % Uncorrected results for the diagonal
    result.uncorr.diag.H                       = diag(squeeze(H_uncorr));
    result.uncorr.diag.P                       = diag(squeeze(P_uncorr));
    result.uncorr.diag.CItwosided(1,:,:)       = diag(squeeze(CItwosided(1,:,:)));
    result.uncorr.diag.CItwosided(2,:,:)       = diag(squeeze(CItwosided(2,:,:)));
    result.uncorr.diag.CItwosided_diff         = diff(result.uncorr.diag.CItwosided);
    result.uncorr.diag.CItwosided_diff_max     = max(diff(result.uncorr.diag.CItwosided));
    result.uncorr.diag.CItwosided_diff_min     = min(diff(result.uncorr.diag.CItwosided));
    result.uncorr.diag.CItwosided_diff_mean    = mean(diff(result.uncorr.diag.CItwosided));
    
    %% Bonf & holm corrections
    dispv(1, '  B. Bonferroni and Holm-Bonferroni multiple comparison corrections...');
       
    % p-value divided by number of comparisons (45 timepoints)
    bonf_alpha = base_alpha/n_tests/n_tests;
    dispv(1, '    Bonferroni t-test alpha set at %g.', bonf_alpha);
    [H_bonf, bonf_P, ~, ~] = ttest(result.ACC_elements, ...
        chancelvl, 'alpha', bonf_alpha, 'tail', 'right');
    H_bonf(bonf_P>(1-bonf_alpha)) = -1;
    isnull = squeeze(round(bonf_P, 4)) == 1 & (squeeze(round(result.mean, 4)) == ...
        round(chancelvl, 4) == 1); % equal to chance
    H_bonf(isnull == 1) = 0;
    
    % Treat NaNs (likely due to extremely similar or exactly the same
    % values) as zero
    if any(isnan(H_bonf), 'all')
        H_bonf(isnan(H_bonf)) = 0;
    end
    
    % Two-sided for CI
    [~, ~,result.bonf.CItwosided, ~]    = ttest(result.ACC_elements, chancelvl, 'alpha', bonf_alpha); 
    
    % Store
    result.bonf.H                       = squeeze(H_bonf);
    result.bonf.alpha                   = bonf_alpha;
    result.bonf.CItwosided_min_meanval  = result.bonf.CItwosided(:, min_ind);
    result.bonf.CItwosided_max_meanval  = result.bonf.CItwosided(:, max_ind);
    result.bonf.CItwosided_diff         = squeeze(diff(result.bonf.CItwosided));
    result.bonf.CItwosided_diff_max     = squeeze(max(diff(result.bonf.CItwosided)));
    result.bonf.CItwosided_diff_min     = squeeze(min(diff(result.bonf.CItwosided)));
    result.bonf.CItwosided_diff_mean    = squeeze(mean(diff(result.bonf.CItwosided)));
    
    % Corrected results for the diagonal
    result.bonf.diag.H                       = diag(squeeze(H_bonf));
    result.bonf.diag.P                       = diag(squeeze(bonf_P));
    result.bonf.diag.CItwosided(1,:,:)       = diag(squeeze(result.bonf.CItwosided(1,:,:)));
    result.bonf.diag.CItwosided(2,:,:)       = diag(squeeze(result.bonf.CItwosided(2,:,:)));
    result.bonf.diag.CItwosided_diff         = diff(result.bonf.diag.CItwosided);
    result.bonf.diag.CItwosided_diff_max     = max(diff(result.bonf.diag.CItwosided));
    result.bonf.diag.CItwosided_diff_min     = min(diff(result.bonf.diag.CItwosided));
    result.bonf.diag.CItwosided_diff_mean    = mean(diff(result.bonf.diag.CItwosided));
    
    % Holm diagonal
    [result.holm.H_diag, result.holm.nullrejections, ...
        result.holm.info] = compute_bonf_holm_cor(result.uncorr.diag.P, ...
        base_alpha, 1, result.diag.mean, chancelvl);
    
    % Holm column-wise (test) matrix
    result.holm.H = compute_bonf_holm_cor(result.uncorr.P(:), ...
        base_alpha, 1, result.mean(:), chancelvl);
    result.holm.H = reshape(result.holm.H, size(result.uncorr.P));
    
    %% FDR
    dispv(1, '\n  C. FDR correction using B&H algorithm fit for correlated results...\n');
    
    % --------------- %
    result = fdr_matrix_mcc(result, chancelvl, base_alpha, 0.001);
    % --------------- %
   
    %% Sign permutation test (size and t-mass)
    dispv(1, '\n  D. Sign permutation tests...\n');
    
    result.signperm.param = [];
    result.signperm.param.do_lowertail = true;
    result.signperm.param.max_perm_p = result.uncorr.alpha; 
    result.signperm.param.param.display_each_sec = 30;
    
    % Data quadrants (sign perm test sometimes fails on negative clusters
    % otherwise)
    data_min_chance = shiftdim(result.ACC_elements, 1)-result.chancelevel;
    task_end = 27; rest_start = task_end + 1; % Timing
    task_HRF_cover = task_end + 6;
    
    dq(1).data = data_min_chance;
    dq(1).name = 'full-data';
    
    dq(2).data = data_min_chance(1:task_HRF_cover, 1:task_HRF_cover, :);
    dq(2).name = 'task-task';
    
    dq(3).data = data_min_chance(1:task_end, rest_start:end, :);
    dq(3).name = 'task-rest';
    
    dq(4).data = data_min_chance(rest_start:end, 1:task_end, :);
    dq(4).name = 'rest-task';
    
    dq(5).data = data_min_chance(rest_start:end, rest_start:end, :);
    dq(5).name = 'rest-rest';
    
    for q_ind = 1:length(dq)
        dispv(1, '\n    Sign permutation for "%s" part of the accuracy matrix\n', dq(q_ind).name);
        
        % Do Kai's sign permutation test
        [result.signperm.dq(q_ind).param, ~, ... % result.signperm.p_elementwise
        result.signperm.dq(q_ind).p_cluster_all, ...
        result.signperm.dq(q_ind).cluster_id_map,...
        result.signperm.dq(q_ind).allclusterstat, ~, ...
        result.signperm.dq(q_ind).p_cluster_all_neg, ...
        result.signperm.dq(q_ind).cluster_id_map_neg, ...
        result.signperm.dq(q_ind).allclusterstat_neg, ~, ...
        result.signperm.dq(q_ind).groupmean, ...
        result.signperm.dq(q_ind).groupstats, ...
        ~, ~, ...
        ~, ~, ...
        ~, ~, ...
        result.signperm.dq(q_ind).p_clustermass_all, ...
        result.signperm.dq(q_ind).p_clustermass_all_neg] = signperm_test(dq(q_ind).data, result.signperm.param);

        % Make sig. H map
        % size
        H = double(result.signperm.dq(q_ind).p_cluster_all < result.signperm.param.max_perm_p);
        H_below = result.signperm.dq(q_ind).p_cluster_all_neg > 1-result.signperm.param.max_perm_p;
        H(H_below) = -1; 
        result.signperm.dq(q_ind).H_size = H; % the original permutation test with size only
        result.signperm.dq(q_ind).alpha = result.signperm.param.max_perm_p;
                
        % mass
        H = double(result.signperm.dq(q_ind).p_clustermass_all < result.signperm.param.max_perm_p);
        H_below = result.signperm.dq(q_ind).p_clustermass_all_neg > 1-result.signperm.param.max_perm_p;
        
        H(H_below) = -1; 
        result.signperm.dq(q_ind).H_mass = H; % the probably more sensitive permutation test with sum(t) (but probably still needs the quandrants, because it's simply different task phases)
        result.signperm.dq(q_ind).alpha = result.signperm.param.max_perm_p;
        
        % Store data and name
        result.signperm.dq(q_ind).data = dq(q_ind).data;
        result.signperm.dq(q_ind).name = dq(q_ind).name;
    end
    
    % Reconstruct based on data quadrants
    htype = {'H_mass', 'mass'; 
        'H_size', 'size'};
    for h_ind = 1:size(htype, 1)
        curr_H = htype{h_ind, 1};
        curr_name = htype{h_ind, 2};
        
        if ~isfield(result.signperm.dq(2), curr_H)
            disp([' Skipping reconstructing quandrants for ' curr_name ' because field result.signperm.dq(2).' curr_H ' does not exist (probably not computed).'])
            continue
        end
        sz = size(data_min_chance);
        result.signperm.dq(6).(curr_H) = zeros(sz(1:2)); % pre-allocate
        if task_HRF_cover == task_end % no extended time block
            result.signperm.dq(6).(curr_H)(1:task_end, 1:task_end) = result.signperm.dq(2).(curr_H);
            result.signperm.dq(6).(curr_H)(1:task_end, rest_start:end) = result.signperm.dq(3).(curr_H);
            result.signperm.dq(6).(curr_H)(rest_start:end, 1:task_end) = result.signperm.dq(4).(curr_H);
            result.signperm.dq(6).(curr_H)(rest_start:end, rest_start:end) = result.signperm.dq(5).(curr_H);
        else
            result.signperm.dq(6).(curr_H)(1:task_HRF_cover, 1:task_HRF_cover) = result.signperm.dq(2).(curr_H);
            H_neg = zeros(sz(1:2)); % pre-allocate
            H_neg(1:task_end, rest_start:end) = result.signperm.dq(3).(curr_H);
            H_neg(rest_start:end, 1:task_end) = result.signperm.dq(4).(curr_H);
            result.signperm.dq(6).(curr_H) = result.signperm.dq(6).(curr_H) + H_neg; % Warning: will cancel out above/below chance on same element (-1 + 1 = 0)
            result.signperm.dq(6).(curr_H)(rest_start:end, rest_start:end) = result.signperm.dq(5).(curr_H);
        end
        result.signperm.dq(6).name = ['full matrix reconstructed from data quadrant tests (mass) and/or (size)'];
        
        % Reconstruct based on original + data quadrants
        result.signperm.dq(7).(curr_H) = result.signperm.dq(1).(curr_H);
        H_neg = zeros(sz(1:2)); % pre-allocate
        H_neg(1:task_end, rest_start:end) = result.signperm.dq(3).(curr_H);
        H_neg(rest_start:end, 1:task_end) = result.signperm.dq(4).(curr_H);
        result.signperm.dq(7).(curr_H) = result.signperm.dq(7).(curr_H) + H_neg; % Warning: will cancel out above/below chance on same element (-1 + 1 = 0)
        result.signperm.dq(7).name = ['full matrix reconstructed from original result matrix + data quadrant tests (mass) and/or (size)'];
    end
    %% TFCE
    dispv(1, '\n  D(2). TFCE...\n');

    try % open one parpool for all TFCE runs
        warning('off', 'matlab_tfce:workers_discovered') % kais changed version, else has no effec
        parpool(20); % auto shutdown after some minutes
    end
    for q_ind = 1:length(dq)
        dispv(1, '\n    TFCE for "%s" part of the accuracy matrix\n', dq(q_ind).name);
        
        % add singleton z dimension as third dimension
        imgs = [];
        for s_ind = 1:size(dq(q_ind).data, 3)
            imgs(:, :, 1, s_ind) = dq(q_ind).data(:, :, s_ind);
        end
        
        [pcorr_pos2,pcorr_neg2] = matlab_tfce('onesample',2,imgs);
        
        result.tfce.dq(q_ind).pcorr_pos2 = pcorr_pos2;
        result.tfce.dq(q_ind).pcorr_neg2 = pcorr_neg2;
        result.tfce.dq(q_ind).alpha = .05;

        % Make sig. H map
        % size
        H = double(result.tfce.dq(q_ind).pcorr_pos2 < result.tfce.dq(q_ind).alpha);
        H_below = result.tfce.dq(q_ind).pcorr_neg2 > 1-result.tfce.dq(q_ind).alpha ...
            & result.tfce.dq(q_ind).pcorr_neg2 ~= 1; % 1.0 stands for: not tested
        H(H_below) = -1; 
        result.tfce.dq(q_ind).H = H;

        % Store data and name
        result.tfce.dq(q_ind).data = dq(q_ind).data;
        result.tfce.dq(q_ind).name = dq(q_ind).name;
    end
    
    % copied from above
    % Reconstruct based on data quadrants
    sz = size(data_min_chance);
    result.tfce.dq(6).H = zeros(sz(1:2)); % pre-allocate
    if task_HRF_cover == task_end % no extended time block
        result.tfce.dq(6).H(1:task_end, 1:task_end) = result.tfce.dq(2).H;
        result.tfce.dq(6).H(1:task_end, rest_start:end) = result.tfce.dq(3).H;
        result.tfce.dq(6).H(rest_start:end, 1:task_end) = result.tfce.dq(4).H;
        result.tfce.dq(6).H(rest_start:end, rest_start:end) = result.tfce.dq(5).H;
        result.tfce.dq(6).name = 'full matrix reconstructed from data quadrant tests (TFCE)';
    else
        result.tfce.dq(6).H(1:task_HRF_cover, 1:task_HRF_cover) = result.tfce.dq(2).H;
        H_neg = zeros(sz(1:2)); % pre-allocate
        H_neg(1:task_end, rest_start:end) = result.tfce.dq(3).H;
        H_neg(rest_start:end, 1:task_end) = result.tfce.dq(4).H;
        result.tfce.dq(6).H = result.tfce.dq(6).H + H_neg; % Warning: will cancel out above/below chance on same element (-1 + 1 = 0)
        result.tfce.dq(6).H(rest_start:end, rest_start:end) = result.tfce.dq(5).H;
        result.tfce.dq(6).name = 'full matrix reconstructed from HRF-extended data quadrant tests (TFCE)';
    end
    
    % Reconstruct based on original + data quadrants
    result.tfce.dq(7).H = result.tfce.dq(1).H;
    H_neg = zeros(sz(1:2)); % pre-allocate
    H_neg(1:task_end, rest_start:end) = result.tfce.dq(3).H;
    H_neg(rest_start:end, 1:task_end) = result.tfce.dq(4).H;
    result.tfce.dq(7).H = result.tfce.dq(7).H + H_neg; % Warning: will cancel out above/below chance on same element (-1 + 1 = 0)
    result.tfce.dq(7).name = 'full matrix reconstructed from original result matrix + data quadrant tests (TFCE)';
   
    %% Store
    comp_results{a_ind} = result;
    
end % comparison

%% Statistics over all comparisons

% Check that analysis will have >1 comparison
if length(all_mask_ACCs.ACCs) > 1 && length(comp_results) > 1
    %% Pair-wise t-tests to check differences between comparisons
    % anova is unneccesarily complex and will prob violate independence assumptions
    
    % Set timing
    tptask = 2;
    taskend = 27;
    reststart = 28;
    tprow = size(result.ACC_elements, 2);
    
    % Prepare
    mask_total = struct();
    mask_total.mask_name = all_mask_ACCs.mask_name;
    for p_ind = 1:length(comp_results)
        mask_total.all_mask_ACCs(:,:,:,p_ind) = shiftdim(comp_results{p_ind}.ACC_elements, 1);
    end
    
    dispv(1, '  D. Pair-wise t-test: averaged across block phases, only works for three comparisons!');
    for p_ind = 1:length(comp_results)
        
        % Wrap around
        wrap_ind = p_ind + 1; if p_ind == length(comp_results), wrap_ind = 1; end
        
        tt1 = squeeze(mean(mask_total.all_mask_ACCs(tptask:taskend, tptask:taskend, :,p_ind), [1 2]));
        tr1 = squeeze(mean(mask_total.all_mask_ACCs(tptask:taskend, reststart:tprow, :,p_ind), [1 2]));
        rt1 = squeeze(mean(mask_total.all_mask_ACCs(reststart:tprow, tptask:taskend, :,p_ind), [1 2]));
        rr1 = squeeze(mean(mask_total.all_mask_ACCs(reststart:tprow, reststart:tprow, :,p_ind), [1 2]));
        tt2 = squeeze(mean(mask_total.all_mask_ACCs(tptask:taskend, tptask:taskend, :,wrap_ind), [1 2]));
        tr2 = squeeze(mean(mask_total.all_mask_ACCs(tptask:taskend, reststart:tprow, :,wrap_ind), [1 2]));
        rt2 = squeeze(mean(mask_total.all_mask_ACCs(reststart:tprow, tptask:taskend, :,wrap_ind), [1 2]));
        rr2 = squeeze(mean(mask_total.all_mask_ACCs(reststart:tprow, reststart:tprow, :,wrap_ind), [1 2]));
        value_set1 = {tt1, tr1, rt1, rr1};
        value_set2 = {tt2, tr2, rt2, rr2};
        
        pairnames{p_ind} = sprintf('%s against %s', comp_results{p_ind}.comparison_str, ...
            comp_results{wrap_ind}.comparison_str); %#ok<*AGROW>
        dispv(1, '    Current pair: %s', pairnames{p_ind});
        
        % Compute t-test with bonf for pairs
        for v_ind = 1:length(value_set1)
            [pair_H(v_ind, p_ind, :), pair_P(v_ind, p_ind, :), ...
                pair_CI(v_ind, p_ind, :), STATS_pair(v_ind, p_ind, :)] = ...
                ttest(value_set1{v_ind}, value_set2{v_ind}, ...
                'alpha', 0.05/(length(value_set1) + length(comp_results)));
            curr_SD     = STATS_pair(v_ind, p_ind).sd;
            curr_size   = sqrt(size(reshape(comp_results{p_ind}.mean.', [], 1), 1)); % 45
            mask_total.pair_stderr(v_ind, p_ind) = curr_SD/curr_size;
        end
        
        mask_total.taskpair_valuediff(:, p_ind) = tt1 - tt2;
        mask_total.taskpair_meandiff(:, p_ind)  = mean(mask_total.taskpair_valuediff(:, p_ind));
        mask_total.valuenames = {'Task-task decoding', 'Task-rest decoding', ...
            'Rest-task decoding', 'Rest-rest decoding'};
    end
    
    mask_total.STATS_pair   = STATS_pair;
    mask_total.pair_H       = pair_H;
    mask_total.pair_P       = pair_P;
    mask_total.pair_CI      = pair_CI;
    mask_total.pairnames    = pairnames;
else
    % Place single comparison as mask_total
    mask_total              = comp_results{1};
    mask_total.n_comp       = 1;
end

%% Get binomial randomized distribution and polyfit for current data

% % Check flag for comparisons
% if (isfield(mask_total, 'n_comp') && mask_total.n_comp > 1) || length(comp_results) > 1
%     
%     % Shift subject last
%     maskdata = shiftdim(mask_total.ACC_subs, 1);
%     
% elseif length(comp_results) == 1
%     
%     % Check dims
%     size_check  = size(mask_total.ACC_elements);
%     size_which  = find(size_check == length(subs_todo));
%     if size_which == 1
%         maskdata = shiftdim(mask_total.ACC_elements, 1);
%     else
%         maskdata = mask_total.ACC_elements;
%     end
% end
% 
% % Generate some binomial drawn data
% try
%     [mean_rsquare, goodness, params] = generate_binom_curve(maskdata, ...
%         chancelvl, 1);
% catch
%     dispv(1, '    Binomial fitting failed, possibly due to different data dimensions');
% end
% 
% % Mean stats
% mask_total.binom.mean_rsquare   = mean_rsquare;
% mask_total.binom.goodness       = goodness;
% mask_total.binom.params         = params;

end % func


%% Subfunctions


%% Get measure of binomial fit of distribution given null/mean data
%
% [mean_rsquare, goodness, params] = generate_binom_curve(inputdata, ...
%   chancelevel, distfit)
%
% Will generate binomially distributed, timeresolved train x test data with
% the same size as inputdata (n*n*x matrix).
%
% Here, inputdata is only expected to be accuracy data (not
% accuracy-minus-chance).
%
% Will also attempt curve fitting. P-hat will be the p*1 estimate of n * n
% trials.
%
% Ingmar 14-09-20

function [mean_rsquare, goodness, params] = generate_binom_curve(inputdata, chancelevel, distfit)

dispv(1, '    E. Generating binomial data for curve fitting...');

% Checks
inputsize = size(inputdata);
assert(inputsize(1) == inputsize(2), 'First two dimensions should be equal');
trials = inputsize(1) * inputsize(2);

if ~exist('chancelevel', 'var') || isempty(chancelevel)
    chancelevel = 0.05;
elseif chancelevel > 1 % likely a percentage
    chancelevel = chancelevel/100;
    dispv(1, '    Chancelevel is >1, therefore set as percentage');
end

% Generate binomial data
binom_rand_dist = (binornd(trials, chancelevel, inputsize)./trials) .* 100;

% Fit distribution to mean
if ~exist('distfit', 'var') || isempty(distfit)
    distfit = 0;
elseif distfit > 1
    error('Input variable "distfit" should be a logical');
elseif distfit
    binom_demeaned  = binom_rand_dist - chancelevel * 100;
    binom_rand_dist = binom_demeaned .* std(inputdata, 1, 3) + mean(inputdata, 3);
end

% Get means
mean_input = squeeze(mean(inputdata, 3));
mean_binom = squeeze(mean(binom_rand_dist, 3));

% Do fitting
for r_ind = 1:size(mean_input, 1)
    [~, goodness(r_ind), params(r_ind)] = fit(mean_input(r_ind, :)', ...
        mean_binom(r_ind, :)', 'poly2');
end
mean_rsquare = mean([goodness.rsquare]); % Get R2

dispv(1, '    2nd-degree polyfit to random binomial curve: mean R^2 = %.4f', ...
    mean_rsquare);

end % func


%% NOT USED

%% MASK TOTAL STATS
%     % Set
%     mask_total.n_comp       = length(comp_results);
%     mask_total.ACC_elements = shiftdim(curr_ACC_elements, 1); % comp as 3rd dim
%     mask_total.ACC_subs     = shiftdim(squeeze(mean(mask_total.ACC_elements, 3)), 2); %subs first
%
%     % Regular stats
%     mask_total.mean_ACCs       = squeeze(mean(mask_total.ACC_elements, [3, 4])); % [3 4] = Comparisons and subs
%     mask_total.mean_diag       = diag(mask_total.mean_ACCs);
%     mask_total.sd_all          = mean(std(mask_total.ACC_elements, 0, [3, 4]), 'all');
%     mask_total.var_all         = mean(var(mask_total.ACC_elements, 0, [3, 4]), 'all');
%     mask_total.kurtosis_all    = mean(kurtosis(mask_total.ACC_elements, 1, [3 4]), 'all');
%     mask_total.skewness_all    = mean(skewness(mask_total.ACC_elements, 1, [3 4]), 'all');
%
%     % Get measure of data distribution of full ACC matrix (across subs)
%     mask_total.sd_ACC_elements      = std(mask_total.ACC_elements, 0, [3 4]);
%     mask_total.var_ACC_elements     = var(mask_total.ACC_elements, 0, [3 4]);
%     mask_total.kurtosis_elements    = kurtosis(mask_total.ACC_subs, 1, 3);
%     mask_total.skewness_elements    = skewness(mask_total.ACC_subs, 1, 3);
%     mask_total.kurtosis_minmax      = [min(min(mask_total.kurtosis_elements)), ...
%         max(max(mask_total.kurtosis_elements))];
%     mask_total.skewness_minmax      = [min(min(mask_total.skewness_elements)), ...
%         max(max(mask_total.skewness_elements))];
%
%     % get Bonferroni-corrected stats across comparisons
%     [mask_bonf_H, mask_bonf_P, ~, ~] = ttest(mask_total.ACC_subs, ...
%         chancelvl, 'alpha', bonf_alpha, 'tail', tail);
%     [~, ~,bonf_CItwosided, ~] = ttest(mask_total.ACC_subs, ...
%         chancelvl, 'alpha', bonf_alpha); % two-sided for CI
%     mask_bonf_H(mask_bonf_P>(1-bonf_alpha)) = -1;
%
%     % Bonf
%     mask_total.bonf.H           = squeeze(mask_bonf_H);
%     mask_total.bonf.P           = squeeze(mask_bonf_P);
%     mask_total.bonf.alpha       = bonf_alpha;
%     mask_total.bonf.CItwosided  = bonf_CItwosided;
%     mask_total.bonf.CItwosided_diff = diff(bonf_CItwosided);
%     mask_total.bonf.CItwosided_mean = mean(mask_total.bonf.CItwosided_diff);
%
%     % Holm diagonal
%     [mask_total.holm.H_diag, mask_total.holm.null_reject, mask_total.holm.info] = ...
%         compute_bonf_holm_cor(diag(mask_total.bonf.P), base_alpha, ...
%         1, mask_total.mean_diag, chancelvl);
%
%     % Holm column-wise (test) matrix
%     for c_ind = 1:n_tests
%         mask_total.holm.H(:, c_ind) = compute_bonf_holm_cor(mask_total.bonf.P(:), ...
%             base_alpha, 1, mask_total.mean_ACCs(:), chancelvl);
%     end

%% Store data in excel and csv
%     % Store result of comparison as struct (later put in cell array, same
%     % data structure as all_ACCs)
%     comp_results{a_ind} = result;
%
%     % Store for results table
%     result_to_write.chancelevel                 = result.chancelevel;
%     result_to_write.tail                        = result.tail;
%     result_to_write.min                         = result.min;
%     result_to_write.max                         = result.max;
%     result_to_write.kurtosis                    = result.kurtosis_all;
%     result_to_write.skewness                    = result.skewness_all;
%     result_to_write.diag_min                    = result.diag.min;
%     result_to_write.diag_max                    = result.diag.max;
%     result_to_write.diag_Tminval                = result.diag.Tminval;
%     result_to_write.diag_Tmaxval                = result.diag.Tmaxval;
%     result_to_write.diag_stderr_minval          = result.diag.stderr_minval;
%     result_to_write.diag_stderr_maxval          = result.diag.stderr_maxval;
%     result_to_write.diag_bonf_CItwosided_min_meanval = result.bonf.diag.CItwosided_min_meanval;
%     result_to_write.diag_bonf_CItwosided_max_meanval = result.bonf.diag.CItwosided_max_meanval;
%     result_to_write.diag_bonf_CItwosided_diff_max    = result.bonf.diag.CItwosided_diff_max;
%     result_to_write.diag_bonf_CItwosided_diff_min    = result.bonf.diag.CItwosided_diff_min;
%     result_to_write.diag_bonf_CItwosided_diff_mean   = result.bonf.diag.CItwosided_diff_mean;
%
%     %% Write up comparison result
%
%     struct_table    = struct2table(result_to_write, 'AsArray', 1);
%     file_table      = fullfile(targetfolder, sprintf('%s_comp-%i_result-table%s.csv', ...
%         curr_maskname, a_ind));
%     writetable(struct_table, [file_table '.xlsx'], 'FileType', 'spreadsheet');
%     writetable(struct_table, [file_table '.csv']); % txt

