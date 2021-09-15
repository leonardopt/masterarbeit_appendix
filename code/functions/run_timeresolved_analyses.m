% function run_timeresolved_analyses(cfg, do)
%
% Caller script for related timeresolved analysis functions.
% Will execute full timeresolved decoding/similarity analysis for all given
% subjects and data types. Important input and function calls are
% bordered with % --------- %.
%
% REQUIRES:
%   cfg struct with valid fields, as checked by the
%    check_passed_cfg subfunction (see end)
%   do-struct with the fields used in this function:
%     do.modelgeneration          = 0/1; % simulations
%     do.firstlevel               = 0/1;
%     do.decoding                 = 0/1;
%     do.decoding_grouplevel      = 0/1;
%     do.parfor                   = 0/1; % turn off for debugging
%     do.all_subs                 = 0/1;
%     do.similarityanalysis       = 0/1;
%     do.similarityanalysis_group = 0/1;
%   Paths set for TDT, SAA simulation functions, and SPM(12) in
%    global_all_paths_private.
%
% OUT: Full decoding results of empirical or simulated dataset
%
% Kai, 02/2020, later expanded by Ingmar and Leonardo, 03/2020 - 06-2021

% HISTORY:
% Joint Leos and Ingmars version, not tested yet. Now solely depends on
%  meta_runfiles - KLI 19.02.2021
% Changes to accomodate similarity analyses in current pipeline - Leonardo
%  & Ingmar, 12.08.2020
% Further reorganisation, ready for simulations and sub/group level
%  decoding - Kai & Ingmar, 19.05.2020
% Reorganisation & name change - Kai, 30.04.2020
% Added process control and more superficial & transparent function calling
%  - Ingmar 28/04

function run_timeresolved_analyses(cfg, do)
%% Passed cfg
% Needs "meta-runfiles" to be passed (e.g. for batching or keeping track
% of multiple simulations). Checks relevant fields here.

valid_cfg = check_passed_cfg(cfg); % See subfunc (no manual input required)
if valid_cfg, run_cfg = cfg; clear cfg; end % Pass fields

% Function specifying other paths for stages (manual changes not required)
run_cfg = set_runfile_paths(run_cfg, do);

%% Print some output

run_cfg.verbose = 1; % 2 is too much output (for loops)

dispv(1, '\nAnalysis processes turned on/off:'); % no dispv before adding TDT
disp(do); pause(0.5); % Display process control settings

dispv(1, '\n%i subject(s) selected for analysis:\n', length(run_cfg.subs_todo));
disp(run_cfg.subs_todo); pause(0.5);

dispv(1, '\n--- Data type set to "%s" ---\n', run_cfg.data_type); pause(0.5);
disp(run_cfg);

% Start timing
run_cfg.time_log.creation_date = datestr(now);
run_cfg.time_log.descr = 'stores timing for different stages';
disp(['Starting analysis: ', run_cfg.time_log.creation_date])

%% -------------------------- Start computation -------------------------- %

% Subject-specific parts of the analysis

if do.parfor && (do.modelgeneration || do.firstlevel || do.decoding || do.similarityanalysis)
    warning('full_pipeline:parfor', ['PARFOR turned on for subject looping! ', ...
        'Compatible with decoding_toolbox_v3.999beta']); pause(2);
    
    verbose         = 0; % output will be mismatched anyway
    run_cfg.verbose = verbose;
    
    % Workers: any whole factor of 54 (total subs) is optimal (e.g. 9 or 18)
    if ~isfield(run_cfg, 'parfor_workers')
        run_cfg.parfor_workers = 18;
    end
    
    try % force parallel pool
        pool = parpool(run_cfg.parfor_workers);
    catch
        pool = run_cfg.parfor_workers;
    end
    
    parfor (sub_idx = 1:length(run_cfg.subs_todo), pool)
        % ------------------------------------------------------------%
        subject_specific_computations(do, run_cfg, run_cfg.subs_todo, sub_idx)
        % ------------------------------------------------------------%
    end
else % NORMAL FOR
    for sub_idx = 1:length(run_cfg.subs_todo)
        % ------------------------------------------------------------%
        subject_specific_computations(do, run_cfg, run_cfg.subs_todo, sub_idx)
        % ------------------------------------------------------------%
    end
end

%% Start group-level statistics

if do.decoding_grouplevel
    dispv(1, 'Computing group-level accuracy matrix for %s ...', run_cfg.post_level1); layout_line_break(1);
    % Group levels -------------------------------------------------------%
    afx_timeresolved_decoding_grouplevel(run_cfg); % optional: pass mask number(s)
    % --------------------------------------------------------------------%
    dispv(1, 'Group-level accuracy matrix computed.'); layout_line_break(1);
end

%% Start similarity analysis

if do.similarityanalysis_group
    dispv(1, 'Start group similarity analysis...')
    % --------------------------------------------------------------------%
    afx_similarity_analysis_grouplevel(run_cfg)
    % --------------------------------------------------------------------%
    dispv(1, 'Group similarity analysis - completed.')
end

%% Newly added: reduce filesize of SPM.mat (can be very large) after simulation FIR1

if strcmp(run_cfg.data_type, 'simulated')
    try
        dispv(1,'Attempting to reduce SPM.mat filesizes...')
        reduce_SPM_mats_from_dirs(run_cfg);
    catch e
        dispv(1,'Reducing SPM.mat filesizes did not work: %s', e.message)
    end
end

%% Print some output

dispv(0, '\nFINISHED %s data analysis: %s\n', run_cfg.data_type, run_cfg.post_level1);

end % Main func




%% ------------------- Subject-specific computations ------------------- %%
%
% Subject-specific parts of the analysis

function subject_specific_computations(do, run_cfg, subs_todo, sub_idx)

curr_cfg        = run_cfg; % reset cfg for current subject
curr_cfg.n_subs = length(subs_todo);
curr_sub        = subs_todo(sub_idx);
curr_cfg.subj   = curr_sub;

%% Step 1: DATA GENERATION
% Get additional data generation settings (for simulation only)

if do.modelgeneration && strcmp(curr_cfg.data_type, 'simulated')
            
    % Make preprocessing (sim) folder if not existent
    preproc_dir = curr_cfg.preproc_dir;
    if ~isfolder(preproc_dir)
        mkdir(preproc_dir); dispv(1, '    Folder "%s" created.\n', preproc_dir);
    end
    % Add sub-xx folder (could this be done in one go?)
    preproc_sub = fullfile(preproc_dir, sprintf('sub-%02i', curr_sub));
    if ~isfolder(preproc_sub)
        mkdir(preproc_sub); dispv(1, '    Subject folder %s created.\n', preproc_sub)
    end
    
    dispv(1, 'Starting: Data generation for subject %0i...', curr_sub); layout_line_break(1);
    % Create data, if data has not been loaded -----------------------%
    [curr_cfg, ~, ~] = generate_single_subject_data(curr_cfg, curr_sub); % ~ output mask names
    % ----------------------------------------------------------------%
    dispv(1, 'Finished: Generating new data for subject %0i.', curr_sub); layout_line_break(1);
    
elseif do.modelgeneration && ~strcmp(curr_cfg.data_type, 'simulated')
    warning('Data generation does not work for anything else than curr_cfg.data_type = ''simulated''. Skipping...')
end

%% Step 2: FIRST LEVEL ()without masks

if do.firstlevel
    if strcmp(run_cfg.data_type, 'simulated') % Currently only made for simulated voxel data
        % First levels FIR setup -------------------------------------%
        compute_level1_FIR_test(curr_cfg, curr_sub);
        % ------------------------------------------------------------%
        dispv(1, '\nFirst-level FIR model computed for subject %0i', curr_sub)
    else
        warning('Empirical/ROI/wholebrain level 1 FIR creation not implemented here but is probably already available. Skipping...')
    end
end

%% Step 3: DECODING

if do.decoding
    dispv(1, '\nStarting timeresolved decoding for subject %0i...', curr_sub); layout_line_break(1);
    % Possible output: [res_fnames, mask_data, cfg, results]----------%
    afx_tdt_subject_timeresolved_decoding(curr_cfg);
    % ----------------------------------------------------------------%
    dispv(1, '\nFinished timeresolved decoding for subject %0i.', curr_sub); layout_line_break(1);
end

%% Step 4: SIMILARITY ANALYSIS (CORRELATIONS)

if do.similarityanalysis
    dispv(1, 'Starting similarity analysis'); layout_line_break(1);
    
    if strcmp(curr_cfg.data_type, 'simulated')
        
        % Do similarity analysis for each layer
        for layer = curr_cfg.saa.data_gen.layer_number
            dispv(1, 'Similarity analysis for participant %i - layer %i\n', ...
                curr_sub, layer)
            % ----------------------------------------------------------------%
            similarity_res{layer} = afx_compute_subject_means(curr_cfg, layer);
            % ----------------------------------------------------------------%
            dispv(1, 'Similarity analysis for participant %i - layer %i - completed \n', ...
                curr_sub, layer)
        end
        
        % Store and save
        for layer = 1:length(similarity_res)
            similarity_cfg(layer) = curr_cfg; %#ok<*AGROW>
            similarity_cfg(layer).saa.similarity_analysis_results = similarity_res{layer};
        end
        save_fname = fullfile(curr_cfg.corr_dir, sprintf('sub-%02i', curr_sub), 'similarity_cfg');
        save(save_fname, 'similarity_cfg');
        
    elseif (strcmp(curr_cfg.data_type, 'wholebrain') || strcmp(curr_cfg.data_type, 'ROI'))
        
        % Check
        if isfield(curr_cfg.design, 'make_pairwise') && curr_cfg.design.make_pairwise == 0
            warning('Though curr_cfg.design.make_pairwise == 0, similarity analysis currently only supports pair-wise analyses. Skipping...');
            return;
        end
        
        % Wholebrain/ROI empirical similarity analysis
        for mask = 1:numel(curr_cfg.all_corr_dir_ROI)
            curr_cfg.corr_dir_ROI = curr_cfg.all_corr_dir_ROI{mask}; % current ROI mask
            
            dispv(1, 'Similarity analysis for participant %i and mask %0i', curr_sub, mask)
            % ----------------------------------------------------------------%
            similarity_res = afx_compute_subject_means(curr_cfg);
            % ----------------------------------------------------------------%
            dispv(1, 'Similarity analysis for participant %i and mask %0i - completed', curr_sub, mask)
            
            % Store and save
            similarity_cfg = curr_cfg;
            similarity_cfg.saa.similarity_analysis_results = similarity_res;
            if strcmp(curr_cfg.data_type, 'wholebrain')
                save_fname = fullfile(curr_cfg.corr_dir, sprintf('sub-%02i', curr_sub), 'similarity_cfg');
            elseif strcmp(curr_cfg.data_type, 'ROI')
                save_fname = fullfile(curr_cfg.corr_dir_ROI, sprintf('sub-%02i', curr_sub), 'similarity_cfg');
            end
            save(save_fname, 'similarity_cfg');
        end
    else
        error('Data type unknown')
    end
    dispv(1, 'Similarity analysis (time-resolved Fisher-Z correlations) for subject %i completed\nResults saved in folder', curr_sub); layout_line_break(2);
end

end



%% Basic checks for a passed cfg (e.g. from batch/meta runfile)
%
% Ingmar, 11-Aug-2020

function valid_cfg = check_passed_cfg(cfg)

fields_to_check = {'subs_todo', 'data_type', 'level1', 'post_level1', ...
    'analysis_type', 'derivative_dir'}; % can add cfg.fields to this list

% Checks basic cfg
for f_ind = 1:length(fields_to_check)
    curr_var = genvarname(['cfg.' fields_to_check{f_ind}]);
    assert(isfield(cfg, fields_to_check{f_ind}) && ~isempty(curr_var), ...
        'Passed cfg invalid: cfg.%s is missing or empty. Use set_all_paths() and specify all relevant fields in cfg.', ...
        fields_to_check{f_ind});
end

assert(isfield(cfg, 'results') && isfield(cfg.results, 'output') && ...
    ~isempty(cfg.results.output), ...
    'Passed cfg invalid: cfg.results.output is missing or empty.');

% Check passed cfg for simulation batch
if strcmpi(cfg.data_type, 'simulated')
    assert(isfield(cfg, 'saa') && ~isempty(cfg.saa), ...
        'cfg.saa was not passed')
    assert(isfield(cfg, 'sim_name') && ~isempty(cfg.sim_name), ...
        'cfg.sim_name was not passed');
    if ~isfield(cfg.saa, 'data_gen')
        warning('A cfg was passed to runfile, but it contained no data_gen parameters');
    end
    if ~isfield(cfg, 'flags') || isempty(cfg.flags)
        warning('A cfg was passed to runfile, but it contained no (overwriting) flags');
    end
end

% Check passed cfg for ROI analysis
if strcmpi(cfg.data_type, 'ROI')
    if ~isfield(cfg, 'ROI_folder') || isempty(cfg.ROI_folder)
        warning('cfg.ROI_folder is empty, check if the cfg input is correct!')
        assert(isfield(cfg.files, 'mask') && ~isempty(cfg.files.mask), ...
            'cfg.ROI_folder and cfg.files.mask were both empty. Check input.');
    end
    if ~isfield(cfg, 'ROI_masks') || isempty(cfg.ROI_masks)
        warning([ 'cfg.ROI_masks is empty, check if the cfg input is correct ', ...
            'or that mask filenames are passed directly to cfg.files.mask'])
        assert(isfield(cfg.files, 'mask') && ~isempty(cfg.files.mask), ...
            'cfg.ROI_masks and cfg.files.mask were both empty. Check input.');
    end
end
valid_cfg = 1; % Passed checks
end

%% function run_cfg = set_runfile_paths(run_cfg)
%
% Explains folder structure associated with AFX empirical and simulation
% decoding.
%
% Kai, May 2020

function run_cfg = set_runfile_paths(run_cfg, do)
%% Append sim name to path names

if ~isfield(run_cfg, 'sim_name') || isempty(run_cfg.sim_name)
    sim_attach = '';
else
    sim_attach = ['_' run_cfg.sim_name];
end

%% Preprocessed data (modelled/fmri) folder, with /sub-99 subfolders

switch run_cfg.data_type
    case 'simulated'
        
        layout_line_break(1); disp('<strong>      Using SIMULATED DATA as input</strong>'); layout_line_break(1);
        
        % If the simulation uses the SPM.mat of real data to create data,
        %   two folders need to be specified here:
        % 1. preproc: where the simulated preproc folder should be
        % 2. a folder where a SPM.mat from real data is (this
        % simulated data: sim-name dir in derivatives
        % analysis/share/corinna_kai/afx/derivatives/sim-timeresolved
        run_cfg.preproc_dir= fullfile(run_cfg.derivative_dir, run_cfg.sim_name);
        if isfolder(run_cfg.preproc_dir)
            warning('full_pipeline:sim_dir_exist', ...
                'Current directory set for simulation already exists, be careful not to accidentally overwrite existing simulations!');
            pause(2);
        end
        
        % Special folder, only for data creation: data containing SPM.mat from
        % real fMRI data (used for model setup)
        % /analysis/share/corinna_kai/afx/derivatives/level1-1/.../func
        run_cfg.level1_fmrispmmat_dir = fullfile(run_cfg.all_paths.data_base_dir, 'derivatives', 'level1-1');
        
        % Simulated similarity analysis
        if do.similarityanalysis || do.similarityanalysis_group
            disp('    Note: correlations are stored slightly differently than decoding.');
            run_cfg.corr_dir = fullfile(run_cfg.derivative_dir, ...
                ['similarity-analysis_from' strrep(run_cfg.level1, 'level1', '') sim_attach]);
        end
        
    case {'wholebrain', 'ROI'}
        
        % Real data: bidsdirectory, e.g. analysis/share/corinna_kai/afx
        % (original directly in BIDS folder, not in derivatives!)
        layout_line_break(1); disp('<strong>      Using REAL DATA as input</strong>');
        layout_line_break(1); pause(2)
        run_cfg.preproc_dir= fullfile(run_cfg.derivative_dir, '');
        
        % Similarity analysis folders
        if strcmp(run_cfg.data_type, 'ROI') && (do.similarityanalysis || do.similarityanalysis_group)
            dispv(1, '    Note: correlations are stored slightly differently than decoding.');
            run_cfg.corr_dir = fullfile(run_cfg.derivative_dir, ...
                ['similarity-analysis_from' strrep(run_cfg.level1, 'level1', '') '_ROI']);
            
            % Loop over ROI masks if multiple are passed
            if ~isfield(run_cfg, 'ROI') && isfield(run_cfg, 'ROI_masks')
                for mask = 1:length(run_cfg.ROI_masks)
                    ROI_name = run_cfg.ROI_masks{mask};
                    
                    % Compatibility
                    if endsWith(ROI_name, '.nii')
                        ROI_name = ROI_name(1:end-4);
                    elseif ~isempty(fileparts(ROI_name))
                        [~, ROI_name, ~] = fileparts(ROI_name);
                    end
                    
                    % Create
                    run_cfg.all_corr_dir_ROI{mask} = fullfile(run_cfg.corr_dir, ROI_name);
                    if ~isfolder(run_cfg.all_corr_dir_ROI{mask}), mkdir(run_cfg.all_corr_dir_ROI{mask}); end
                end
            else
                % Create
                run_cfg.corr_dir_ROI = fullfile(run_cfg.corr_dir, run_cfg.ROI);
                if ~isfolder(run_cfg.corr_dir_ROI), mkdir(run_cfg.corr_dir_ROI); end
            end
        end
        
    otherwise
        error('Unknown datatype!')
end

% Rename data_base_dir as bids_dir
run_cfg.bids_dir = run_cfg.all_paths.data_base_dir;

% run_cfg.level1 directory of current pipeline modelled data / preprocessed fmri data
% e.g. '/analysis/share/corinna_kai/afx/derivatives/level1-FIR1_sim-timeresolved'
run_cfg.level1_dir = fullfile(run_cfg.derivative_dir, [run_cfg.level1 sim_attach]);

% level2 dir with subject specific post-level1 analysis (e.g. decoding or
% correlation) AND group-level results (called level2 in the past)
% new: level2 dir contains from which level1 dir data came + simulation name
% subfolders:
%   /sub-99 (for subject 99)
%   /group  (n: number of subjects in that group level analysis added to result files)
% e.g. afx/derivatives/decoding-timeresolved_from-FIR1_sim-timeresolved'
run_cfg.level2_dir  = fullfile(run_cfg.derivative_dir, ...
    [run_cfg.post_level1  '_from' strrep(run_cfg.level1, 'level1', '')  sim_attach]);

end
