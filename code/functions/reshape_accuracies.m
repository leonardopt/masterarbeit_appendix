%% function reshape_accuracies(cfg, results, curr_comparison, m_ind, maskname)
%
% Transforms data and adds some info to be usable per mask.
%
% Kai and Ingmar, May 2020

% Separated from afx_tdt_subject_timeresolved_decoding.m
% Ingmar 17-03-21
%
% Updated for searchlights with multiple CV designs - Ingmar 07-21

function ACCs = reshape_accuracies(cfg, results, curr_comparison, m_ind, maskname)

% init
ACCs = struct();

% add name of comparison (this can surely be improved)
% get all steps of current comparison
steps_idcs = cfg.design.condition_list_index_per_step == curr_comparison;

%% --- Check which analysis ---

% for  joint designs of two-class cases: use this code
if (isfield(cfg.design, 'make_pairwise') && cfg.design.make_pairwise) || ...
        ~any(strcmp({'standard_cv', 'standard_difficulty_generalization'}, cfg.analysis_type))
    
    % assert that there are only two labels
    assert(cfg.design.n_cond_per_step == 2, 'Only expected 2 labels per decoding step. Check if settings and label input are correct.')
    
    try
        % get current label index
        c1label = unique(cfg.design.conditions_per_step(1, steps_idcs));
        c2label = unique(cfg.design.conditions_per_step(2, steps_idcs));
        
        % assert they are unique
        assert(numel(c1label)==1, ...
            ['Different values for c1val (label of condition 1 in all decoding steps ', ...
            'that should be used for the current ACC combination), should be unique; please check.']);
        assert(numel(c2label)==1, ['Different values for c2val (label of condition 2 in all decoding steps ', ...
            'that should be used for the current ACC combination), should be unique; please check.']);
        
        % get string for condition and save to title
        cond1str = cfg.design.label_val2str_map.common_descr{cfg.design.label_val2str_map.label==c1label};
        cond2str = cfg.design.label_val2str_map.common_descr{cfg.design.label_val2str_map.label==c2label};
        ACCs.title = [cond1str ' vs ' cond2str];
    catch title_err
        warningv('subject_decoding:cond_title_failed', ...
            'Title construction for current decoding comparison (%i) failed: %s.', ...
            curr_comparison, title_err.message)
        ACCs.title = sprintf('Comparison %i', curr_comparison);
    end
    
    % For multiclass (e.g. 3-class) classification: take all labels from training, and test set
elseif any(strcmp({'standard_cv', 'standard_difficulty_generalization'}, cfg.analysis_type))
    
    % Check
    assert(all(all(diff(cfg.design.conditions_per_step, 1, 2)==0)), ...
        'For multiclass, the labels for all decoding stelps should be the same');
    
    ACCs.title = 'All labels: Sim vs. Add vs. Con'; % hardcoded!
    
else
    error('Unexpected combination of parameters. Check settings.');
end

% add maskname
ACCs.maskname  = maskname;

% Add axes labels that do not contain class info, only time info
ACCs.axislabel = cfg.design.conditions2set_map_axislabel{1};

%% Place accuracy data

if strcmp(cfg.data_type, 'searchlight') % ------ SL data ------ %
            
    % Init, store each sl map in 3rd dim
    voxspace_mask   = cfg.searchlight.voxspace_mask;
    ACCs.SL_subset  = cfg.searchlight.subset;
    ACCs.SL_ACC = nan(numel(unique(cfg.design.set2conditions_map(:, 2))), ...
        numel(unique(cfg.design.set2conditions_map(:, 3))), ...
        numel(ACCs.SL_subset(:, 1)));
    ACCs.SL_IND = ACCs.SL_ACC; % same size
    
    % Checks
    assert(isfield(cfg.searchlight, 'voxspace_mask'), ...
        'Provide a mask of logicals indexing the subset of searchlight locations in 3D space in the original volume dimensions.');
    
    if size(ACCs.SL_subset, 1) ~= sum(voxspace_mask, 'all')
        warningv('SL_subset_no_mask_match', ...
            'Searchlight subset points does not match the number of accuracies in output!');
    end
    assert(all(cfg.datainfo.dim == size(voxspace_mask)), ...
        'Specifiek voxel mask does not agree with dimensions in cfg.datainfo!');
    if size(results.accuracy.set(1).output, 1) ~= size(ACCs.SL_subset, 1)
        warningv('cfg.datainfo_no_mask_match', ...    
        'Searchlight accuracies and the chosen voxel subset do not match in length! Please check.');
    end
    
    % n searchlight map values per train x test timepoint combo
    for set_ind = 1:length([results.accuracy.set.set_id])
        
        % get output for current set
        nonzero_sl_vals = nonzeros(results.accuracy.set(set_ind).output);
        
        % cfg.design.set2conditions_map is comparison x train x test point
        % assign multiple variables to cell array in one line
        curr_set = results.accuracy.set(set_ind).set_id; % is this redundant?
        curr_map = num2cell(cfg.design.set2conditions_map(curr_set, :)); 
        [~, train_ind, test_ind] = curr_map{:}; % unpack cell array
        
        % Ordered by linear index in 3rd dim
        for sl_ind = 1:length(nonzero_sl_vals)
            ACCs.SL_ACC(train_ind, test_ind, sl_ind) = nonzero_sl_vals(sl_ind);          
        end
    end
    
else % ------ no SL data ------ %
    
    % Init
    ACCs.ACC = nan(numel(unique(cfg.design.set2conditions_map(:, 2))), ...
        numel(unique(cfg.design.set2conditions_map(:, 3))));

    % 1 value per train x test timepoint combo
    for set_ind = 1:length([results.accuracy.set.set_id])
        
        % get output for current set for current mask
        curr_val = results.accuracy.set(set_ind).output(m_ind);
        curr_set = results.accuracy.set(set_ind).set_id;
        
        % cfg.design.set2conditions_map is comparison x train x test point
        % assign multiple variables to cell array in one line
        curr_map = num2cell(cfg.design.set2conditions_map(curr_set, :)); 
        [comp_ind, train_ind, test_ind] = curr_map{:}; % unpack cell array
        ACCs.ACC(train_ind, test_ind)   = curr_val;
    end
    
    % Check for NaNs
    if any(isnan(ACCs.ACC))
    warning('Nans in ACCs{%i}.ACC: %i nans in %i values', comp_ind, ...
        sum(isnan(ACCs{comp_ind}.ACC(:))), numel(ACCs{comp_ind}.ACC));
    end
end

end
