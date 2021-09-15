% function afx_timeresolved_decoding_grouplevel(run_cfg)
%
% Do group analysis for time resolved decoding analyses from
%   afx_tdt_subject_timeresolved_decoding.m
% called by
%   run_timeresolved_analyses.m
%
% This function uses cell structures to include all results per mask and
% per subsequent comparison, all of which can vary per analysis.
%
% IN
%   run_cfg: from pipeline, containing
%   run_cfg.subs_todo: subjects to include for analysis (n subjs might be
%       added at result folder name if unexpected n)
%   run_cfg.level2_dir: directory that contains decoding results in sub-99
%       folders for each subject. Results will be written to 'group'
%       subdirectory.
%
% OUT
%   group result files per mask in fullfile(run_cfg.level2_dir, 'group')
%   resultfiles.mats: cellstr list of resulting matfile names
%   resultfiles.figures: % cellstr list of resulting figure names
% Kai, 2020-05-19

% Added timeseries/timeslice plotting functions, Ingmar 2020-06-02

function afx_timeresolved_decoding_grouplevel(run_cfg)
%% Control

do_translate    = 1; % translate German names to English
do.clims        = 1; % equal colorbar limits for result plots

%% Get info

subs_todo    = run_cfg.subs_todo;
total_subs   = length(subs_todo);
sourcefolder = run_cfg.level2_dir;
targetfolder = fullfile(sourcefolder, 'group');
if ~isfolder(targetfolder), mkdir(targetfolder); end

if total_subs < 3
    warning('Minimum sample size for a t-test comparison is n=3. Functions will likely give errors.')
elseif total_subs < 7
    warning('grouplevel_decoding:small_sample', ...
        ['Group level decoding results do not mean anything for small ', ...
        'sample sizes, and plotting will likely be off']);
end

if ~exist('masks_ind_todo', 'var')
    dispv(1, ['Group-level decoding: masks_ind_todo missing as input, ', ...
        'will try to retrieve mask data automatically.']);
end

%% First, check if data present in group folder

dir_output = dir(fullfile(targetfolder, 'decoding-groupstatistics*'));

if ~isempty({dir_output.name})
    dir_names = {dir_output.name};
    assertstr = ['More than 1 group level save was found in %s, but the group-level pipeline can only process one at a time. ', ...
        'Choose one group level save for the current analysis and (re)move the other.'];
    assert(length(dir_names) == 1, assertstr, targetfolder); clear assertstr;
    dispv(1, 'Group-level decoding: "%s" found.', dir_output.name);
    
    try
        savefile = fullfile(targetfolder, dir_names{:});
        cacheddata = load(savefile, 'all_results', 'run_cfg');
        
        disp(' ')
        disp('Loading precomputed grouplevel results')
        disp(' ')        
        disp('Warning: No check implemented if the current run_cfg is consistent with the loaded run_cfg')
        disp('If you want to used the cached run_cfg (will overwrite the set run_cfg):')
        disp('  run_cfg = cacheddata.run_cfg;' )
        disp('  dbcont');
        disp(' ')
        disp('Call only dbcont if you want to use the run_cfg that you defined with the loaded data')
        disp(' ')
        disp('If you DONT KNOW what you do or if the current run_cfg fits the loaded run_cfg:')
        disp('  Call:: error()');
        disp('Yes as simple as that - this will recompute the group level analyses')
        keyboard
        
        all_results = cacheddata.all_results;
    catch e
        warning('Loading cached group levels failed, starting new subject-level loading (see errormsg below)...');
        disp(['    Error: ' e.message]);
    end
    dispv(1, ['\nI. Data already loaded, skipping data retrieval. \n    ', ...
        'Delete/move the found groupstatistics .mat file to restart data compilation.']);
    run_cfg.groupresults_present = 1;
    
    dispv(1, '    Stopping run to inspect results (see "all_results" struct). Type dbcont to continue to plotting.');
    keyboard % to jump into current variables
end % .mat present

%% Load data for all subjects

if ~isfield(run_cfg, 'groupresults_present')
    dispv(1, '\nI. Starting data retrieval for all given subjects (n=%0i)\n', total_subs);
    
    % Some settings
    run_cfg.disable_grouploading_checks = 1; % disable checks for now
    run_cfg.do_parfor                   = 1; % parallel loading
    run_cfg.group_template              = 1; % get template
    
    % Get data for first subject as comparison
    [~, ~, first_data, ~] = afx_timeresolved_decoding_load_group_data(run_cfg, 1, sourcefolder);
    dispv(1, '    Data loaded for first subject as template...\n');
    run_cfg.group_template = 0; % real data
    
    % Initialize
    load_SNR    = cell(size(total_subs));
    load_ACC    = cell(size(total_subs));
    load_CONF   = cell(size(total_subs));
    
    % Get for each subject and compare to first subject (include first subject to be sure all works)
    if run_cfg.do_parfor && total_subs > 3
        
        % Set workers
        parfor_workers = 9;
        if total_subs < 24
            try 
                parfor_workers = parpool(total_subs); 
            end
        else
            try 
                parfor_workers = parpool(20);
            end
        end
        
        % Optimize parfor
        run_cfg_parfor      = parallel.pool.Constant(run_cfg);
        first_data_parfor   = parallel.pool.Constant(first_data);
        
        dispv(1,'    Now starting parallel loading...');
        parfor (subj = 1:total_subs, parfor_workers)
            % --------------------------------------------------------------------%
            [load_ACC{subj}, load_CONF{subj}, ~, ~, load_SNR{subj}] = ...
                afx_timeresolved_decoding_load_group_data(run_cfg_parfor.Value, ...
                subj, sourcefolder, first_data_parfor.Value); % put into target structure
            % --------------------------------------------------------------------%
        end
    else
        for subj = 1:total_subs
            % --------------------------------------------------------------------%
            [load_ACC{subj}, load_CONF{subj}, ~, ~, load_SNR{subj}] = ...
                afx_timeresolved_decoding_load_group_data(run_cfg, ...
                subj, sourcefolder, first_data); % put into target structure
            % --------------------------------------------------------------------%
        end
    end
    
    dispv(1, '  B. First-level accuracies (and possible confusion matrices) loaded succesfully.');
    
    %% Assemble target structure for ACCs
    
    all_ACC = cell(numel(load_ACC{1}), 1); % init
    for subj = 1:total_subs
        for mask = 1:numel(load_ACC{1})
            
            % Result vars
            all_results = struct(); % init
            all_results.masks_ind_todo = 1:numel(all_ACC);
            
            % Accuracy matrices, filenames, titles, axes, mask names
            for comp = 1:numel(load_ACC{1}{mask}.ACCs)
                ACCs_todo = load_ACC{subj}{mask};
                all_ACC{mask}.ACCs{comp}.ACC(:, :, subj) = ACCs_todo.ACCs{comp}.ACC(:, :, subj);
                all_ACC{mask}.full_datafile{subj}       = ACCs_todo.full_datafile{subj};
                all_ACC{mask}.ACCs{comp}.title          = ACCs_todo.ACCs{comp}.title;
                all_ACC{mask}.ACCs{comp}.axislabel      = ACCs_todo.ACCs{comp}.axislabel;
                all_ACC{mask}.mask_name                 = ACCs_todo.mask_name;
            end % comparisons
        end % masks
    end % subs
    dispv(1, '  C. ACC data sorted in structs per subject, mask, and comparison.');
    
    %% Assemble confusion matrix data
    
    if exist('load_CONF', 'var') && ~isempty(load_CONF{1}) && ~isempty(load_CONF{end}) % check
        
        % Compatbility
        for mask = 1:numel(load_CONF{1})
            for subj = 1:total_subs
                if numel(load_CONF{subj}{mask}.CONFMAT) == 1 % So many cells
                    load_CONF{subj}{mask}.CONFMAT = load_CONF{subj}{mask}.CONFMAT{:};
                end
            end
        end
        
        % Unpack per sub for curr mask, take mean, place in cell
        for mask = 1:numel(load_CONF{1})
            for el_ind = 1:numel(load_CONF{1}{1}.CONFMAT)
                for subj = 1:total_subs
                    confmat_el(:,:,subj) = load_CONF{subj}{mask}.CONFMAT{el_ind}.MAT;  %#ok<*AGROW>
                    mean_CONFMAT_els{el_ind} = mean(confmat_el, 3);
                end
            end
            CONFMAT{mask} = mean_CONFMAT_els; % store
        end
        
        % Save element indices
        for el_ind = 1:numel(load_CONF{1}{1}.CONFMAT)
            CONFMAT_compstr{el_ind} = load_CONF{1}{1}.CONFMAT{el_ind}.compstr;
            CONFMAT_elstr{el_ind} = load_CONF{1}{1}.CONFMAT{el_ind}.elstr;
        end
        
        % Store in struct
        for mask = 1:numel(load_CONF{1})
            all_results(mask).CONFMAT = CONFMAT{mask};
            all_results(mask).CONFMAT_compstr = CONFMAT_compstr;
            all_results(mask).CONFMAT_elstr = CONFMAT_elstr;
        end
        dispv(1, '  D. CONFMAT data sorted in cells and struct per mask.');
    end % confmat
    
    %% Compute stats
    dispv(1, '\nII. Constructing group-level summary statistics...\n');
    
    warning('Trying parfor, change to for if problems occur')
    length_all_ACC = length(all_ACC);
    
    for mask = 1:length_all_ACC % Loop % ADD parfor for faster computing, not for debugging or devel
        
        % Display
        comp_maskname = all_ACC{mask}.mask_name; % retrieve
        fprintf('Processing grouplevel "%s" (mask %i/%i)\n', comp_maskname, mask, length_all_ACC);
        % Get results and store them ---------------------------------------- %
        [all_results(mask).comp_results, all_results(mask).mask_total] = afx_grouplevel_decoding_stats(run_cfg, ...
            all_ACC{mask}, mask); %#ok<*ASGLU>
        % --------------------------------------------------------------------%
        
        % Check if SNR was calculated (WAS WRONG, removed to allow parfor and we
        % dont use it (kai, 2021-09-06)
    %         if strcmpi(run_cfg.data_type, 'simulated') && exist('loaded_SNR', 'var') ...
    %                 && ~isempty(load_SNR{1}) && ~isempty(load_SNR{end})
    %             
    %             % Assembling SNR vals: relies on exact variable order found in
    %             % load_group_data! done this way as parfor workaround
    %             dispv(1,'    Averaging SNR & CNR of individual results (ignoring cue noise)...');
    %             % implementaton is certainly wrong, because loop over
    %             % total_subs is not used. removing
    %             for subj = 1:total_subs
    %                 for layer = 1:length(load_SNR{1}{1})
    %                     snr_db_total(layer, :)   = load_SNR{layer}{1}; % Model timeseries only %#ok<*AGROW,*NASGU>
    %                     CNR_total(layer, :)      = load_SNR{layer}{2};
    %                 end % layers
    %             end % sub
    %             
    %             % Store
    %             all_results(mask).total_SNR = mean(snr_db_total, 1);
    %             all_results(mask).total_CNR = mean(CNR_total, 1);
    %         end % SNR
    end % mask
else
    % Basic check for confmat
    if isfield(all_results(1), 'CONFMAT')
        assert(~isempty(all_results(1).CONFMAT), 'isempty(all_results(1).CONFMAT == 1')
        assert(~isempty(all_results(end).CONFMAT), 'isempty(all_results(end).CONFMAT == 1')
    end
end % data retrieval and compiling

%% Create result plots

dispv(1, '\nIII. Constructing plots...\n');

% prepare parfor loop
if ~isfield(run_cfg, 'do_MCC'), run_cfg.do_MCC = []; end % use default
if exist('CONFMAT', 'var')
    if_exist_CONFMAT_else_empty = CONFMAT;
else
    if_exist_CONFMAT_else_empty = [];
end

disp('Trying parfor for constructing the plots, change to for debugging or if you experience problems') % DEVEL
parfor mask = 1:length(all_results(1).masks_ind_todo)
    
    % Retrieve
    comp_stats = all_results(mask).comp_results;
    comp_maskname = all_results(mask).mask_total.mask_name;
    
    % Translate
    if do_translate && any(contains(comp_maskname, {'aeh', 'kon', 'sum'}, 'IgnoreCase', true))
        comp_maskname = strrep(comp_maskname, 'aeh', 'sim');
        comp_maskname = strrep(comp_maskname, 'kon', 'con');
        comp_maskname = strrep(comp_maskname, 'sum', 'add');
        comp_maskname = strrep(comp_maskname, 'Aeh', 'Sim');
        comp_maskname = strrep(comp_maskname, 'Kon', 'Con');
        comp_maskname = strrep(comp_maskname, 'Sum', 'Add');
    end
    
    % Make mask name clearer
    if strcmpi(comp_maskname, 'mask')
        comp_maskname = [run_cfg.data_type '_mask']; % e.g. wholebrain_mask
    end
    dispv(1, '    Mask %i: %s', mask, comp_maskname);
    
    for comp = 1:length(comp_stats) % Class comparisons
        
        comp_res = comp_stats{comp}; % Retrieve
        
        % Check for FDR corrected stats
        if ~isfield(comp_res, 'fdr')
            try chancelvl = run_cfg.results.chancelevel;
            catch
                chancelvl = 100/(length(strfind(comp_res.title, 'vs')) + 1); % set
            end
            % -------------- %
            comp_res = fdr_matrix_mcc(comp_res, chancelvl); % Adds FDR-corrected sig. matrix
            % -------------- %
        end
        
        % Specify target folder for saving
        targetfile1 = fullfile(targetfolder, strrep(strrep(strrep(['decoding-ttest_', ...
            comp_res.group_resname], ' ', ''), '.', '_'), ':', ''));
        
        % Extend caxis limits if xclass. Determine the data range across
        % all elements involved by using e.g. the
        % quick_results_script.m in ../ROI_scripts
        if do.clims && contains(run_cfg.analysis_type, {'xclass', 'cross'}) || ...
                strcmpi(run_cfg.data_type, 'ROI')
            % ------------- %
            comp_res.extend_caxis = [20 75; 0 20]; % set manually
            % ------------- %
            dispv(1, '    --- CLIMS --- Manual color limits (comp_result.extend_caxis) for accs. and CIs was set to:');
            dispv(1, comp_res.extend_caxis); layout_line_break(1);
        end
        
        %% Plotting
        
        % Results --------------------------------------------------- %
        dispv(1, '    Constructing decoding result plots');
        [customticks, customlabels] = afx_set_timings(0); % Full axes ticks and labels
        % plot_decoding_results(comp_res, targetfile1, customticks, customlabels, run_cfg.data_type, run_cfg.do_MCC);
        plot_all_MCC_results(comp_res, [targetfile1, '_all-MCC_single-file'], run_cfg.data_type); % 1 big file
        % ----------------------------------------------------------- %
        
        % Confusion matrices for group comp
        dispv(1, '    Creating group confusion matrix for current mask...');
        targetfile2 = fullfile(targetfolder, ['confusion-matrix_n-' num2str(total_subs) '_' comp_maskname]);
        [customticks, customlabels] = afx_set_timings(1); % Fewer axes ticks and labels
        if isfield(all_results(mask), 'CONFMAT') % Compatibility (?)
            curr_CONFMAT = all_results(mask).CONFMAT;
        else
            fprintf('all_results(%i).CONFMAT does not exist. might be normal. Using if_exist_CONFMAT_else_empty\n', mask)
            curr_CONFMAT = if_exist_CONFMAT_else_empty; % for parfor
        end
        if ~isempty(curr_CONFMAT)
            % ----------------------------------------------------------- %
            plotting_type = {'single_matrix'}; % Types: 'subplot' or 'single_matrix'
            plot_group_confusion_matrix(curr_CONFMAT, comp_res, ...
                targetfile2, plotting_type, customticks, customlabels);
            % ----------------------------------------------------------- %
        end
        comp_stats{comp} = comp_res; % Store comparison results back
        close ALL;
    end % comparison
    
    % Store all results back (will be done after all parfor loops have been done)
    all_results_tmp(mask).comp_results = comp_stats;
    all_results_tmp(mask).masks_ind_todo = all_results.masks_ind_todo; % not mask-specific
end % mask
disp('parfor done')

% % Store all results back (will be done after all parfor loops have been done)
for mask = 1:length(all_results(1).masks_ind_todo)
    all_results(mask).comp_results = all_results_tmp(mask).comp_results;
    all_results(mask).masks_ind_todo = all_results_tmp(mask).masks_ind_todo;
end

%% Saving grouplevels

if ~exist('savefile', 'var')
    fname = sprintf('decoding-groupstatistics_n-%02i_%s', total_subs, datestr(now, 'dd-mm-yyyy'));
    save_parfor(run_cfg, all_results, targetfolder, fname); % see subfunc
    dispv(1, '\nIV. Group level results saved as "%s"...\n    Can be found in %s.', fname, targetfolder);
end
dispv(1, 'All grouplevels done & finished plotting.')

end % func



%% Function to enable saving within parfor loop

function save_parfor(run_cfg, all_results, targetfolder, savename)

save_loc_results = fullfile(targetfolder, savename);
try
    save(save_loc_results, 'all_results', 'run_cfg', '-v7.3');
    dispv(1, '    Results saved.');
catch e
    warning('Cannot save group results, see error');
    disp(['    Error: ' e.message]);
end

end

