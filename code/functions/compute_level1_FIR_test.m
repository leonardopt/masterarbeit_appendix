% function compute_level1_FIR_test(cfg, subj)
%
% Sets up first-level analysis (subject-specific) using all possible betas
% (timepoint regressors) to construct a full timeresolved decoding
% analysis. Adapted from older functions such as run_level1_FIR1_test.m and
% its parallel and parfor versions constructed by Kai/Leonardo, ~2019.
%
% IN: the following cfg entries are copied to local variables
%    sourcedir   = cfg.preproc_dir;
%    targetdir   = cfg.level1_dir;
%    tsv_basedir = cfg.bids_dir; % path to subject tsv file, e.g. /analysis/share/corinna_kai/afx/sub-01/func
%    (directly in bids dir, not in derivative) (always the same, not dependent on simulation or not)
%
% in sanity check (if nobody changed that)
%   assert(contains(targetdir, cfg.sim_name), ...)
%
% as flags cfg.flags.overwrite_level1_SPM
%
% if remote setup should be used
%   cfg.use_remote_setup_tmpdir and specific cfg.all_paths dirs (see code
%      below)
%
% subj: subj number as int (should correspond to cfg.subj)
%
% OUT: dir with SPM level 1 data
%  for simulation: *masks*.nii (masking each voxel set),
%   mask_names.mat (with names of these masks)
%   sim_cfg.mat (with cfg of simulation)
%   nii_nam.mat (name of simulation in each voxel)
%     first volume of nii_nam in sim_cfg.nii_model_descr
%     (typically sufficient, because all volumes should have
%     the same content) in sim_cfg.
%
% needs SPM
%
% Kai 2020-05-05

function compute_level1_FIR_test(cfg, subj)
%% Directories and cfg settings

sourcedir   = cfg.preproc_dir;
targetdir   = cfg.level1_dir;
tsv_basedir = cfg.bids_dir; 

dispv(1, '\nStarting compute_level1_FIR_test...\n');

%% Check naming of target dir

dispv(1, '1. Checking that FIR in level1_dir name')
assert(contains(targetdir, 'FIR1'), ...
    'Could not find FIR1 in targetdir, aborting (change if this is a new function)')

if strcmpi(cfg.data_type, 'simulated')
    dispv(1, '2. Checking that simulation name is in target dir')
    assert(contains(targetdir, cfg.sim_name), ...
        'Could not find sim_name %s in targetdir, aborting (change if this is a new function)', cfg.sim_name)
end

%% Current subject sourcedir/targetdir/tsv-design file

sub_str         = sprintf('sub-%02i', subj);
source_dir_subj = fullfile(sourcedir, sub_str);
target_dir_subj = fullfile(targetdir, sub_str);

% check if source directory exists
assert(isfolder(source_dir_subj), 'Source directory %s not found, aborting', source_dir_subj)

% if the source has a /func subdir, use this for source and target
% (currently empirical data has a func subdir)
if exist(fullfile(source_dir_subj, 'func'), 'dir')
    source_dir_subj = fullfile(source_dir_subj, 'func');
    target_dir_subj = fullfile(target_dir_subj, 'func');
end

dispv(1, ['   Sourcedir: ' source_dir_subj])
dispv(1, ['   Targetdir: ' target_dir_subj])

tsv_file = fullfile(tsv_basedir, sub_str, 'func', [sub_str '_task-afx_events_FIR1.tsv']);
% e.g. /analysis/share/corinna_kai/afx/sub-01/func (directly in bids dir, not in derivative)
assert(isfile(tsv_file), 'Could not find design tsv file: %s, please check', tsv_file)
dispv(1, ['   Design tsv file: ' tsv_file]);

%% Start diary

try
    curr_date = strrep(datestr(now), ' ', '_');
    diaryfile = fullfile(target_dir_subj, sprintf('level1_FIR1_output_%s.txt', curr_date));
    diary(diaryfile);
    disp('3. Logging output to diary in subject target folder...')
catch
    warning(['Could not create diary ' diaryfile])
end

%% Skip if data has already been computed
% check if SPM.mat and some betas exist (if so, assume the computation has
% been done)

dispv(1, '4. Checking contents of %s..., ', target_dir_subj);

if isfield(cfg, 'saa') && isfield(cfg.saa.data_gen, 'expected_betacount')
    expected_betas = cfg.saa.data_gen.expected_betacount;
    dispv(1, '   Expected time-resolved beta count: %0i for subject %02i', ...
        expected_betas, subj);
else
    expected_betas = 1850; % value of expected betas depends on amount of runs, ~1900 for 7, ~2150 for 8
    dispv(1, '   No known beta count present in cfg, set at %i.', expected_betas);
end

deletefile_str = '   SPM file deleted, recomputing first levels (FIR model) ';
deletefail_str = 'Deletion failed, will likely be prompted with SPM pop-up: choose to recompute files ';

if exist(fullfile(target_dir_subj, 'SPM.mat'), 'file') == 2
    
    % Get dir list
    beta_dir = dir(fullfile(target_dir_subj, 'beta*.*'));
    if size(beta_dir, 1) >= expected_betas || (size(beta_dir, 1) >= expected_betas - 50)
        dispv(1, '   SPM.mat and (probably) all betas found in %s', target_dir_subj)
        
        % Overwrite
        if isfield(cfg.flags, 'overwrite_level1_SPM') && cfg.flags.overwrite_level1_SPM
            dispv(1, '   RECOMPUTING SPM first level because cfg.flags.overwrite_level1_SPM=true')
            try
                %delete(fullfile(target_dir_subj, 'beta*.*'));
                delete(fullfile(target_dir_subj, 'SPM.mat'));
                dispv(1, deletefile_str);
            catch err
                warning(deletefail_str);
                disp(['    Error: ' err.message]);
            end
        else
            dispv(1, ['   SKIPPING SPM first level because result file already exists. ', ...
                'Set cfg.flags.overwrite_level1_SPM=true to overwrite files'])
            return
        end
    else % likely not enough betas
        dispv(1, '   SPM.mat found, but beta count is lower than expected.')
        if isfield(cfg.flags, 'overwrite_level1_SPM') && cfg.flags.overwrite_level1_SPM
            dispv(1, [deletefile_str, 'since overwrite_level1_SPM == true..., '])
            try
                % delete(fullfile(target_dir_subj, 'beta*.*'));
                delete(fullfile(target_dir_subj, 'SPM.mat'));
                dispv(1, deletefile_str);
            catch err
                warning(deletefail_str);
                disp(['    Error: ' err.message]);
            end
        else % no automatic overwriting
            fname_errorlog = fullfile(targetdir, 'warning_beta_count_subslist.txt');
            write_str = sprintf('%s: Subject %02i failed expected beta count >= %0i. \n', ...
                datestr(now), subj, expected_betas);
            write_str = [write_str, 'Check manually or re-run scripts with overwrite_level1_SPM = 1'];
            try % log in existing file
                fid = fopen(fname_errorlog, 'a+');
                fprintf(fid, ['\n' write_str]); % new line
                fclose(fid);
                dispv(1, '   Logged sub-%02i beta count mismatch to: %s', subj, fname_errorlog);
            catch write_err % probably no writing permissions, though then beta generation will fail anyway
                warning('Could not write errorlog to targetdir:')
                disp(['    Error: ' write_err.message]); 
            end
            
            warning(['Since overwrite_level1_SPM == false, manual beta checking/deletion is required.', ...
                'Output logged in subject diary and across subjects in target folder.'])
            
            % Quit processing as no further runfile iterations will be done
            if isfield(cfg, 'n_subs') && cfg.n_subs == 1 % new
                error('Script is not looping. Aborting.') 
            else
                % Skip, no info available over rest of control processes
                warning('Skipping current subject FIR1 model generation...')
                diary off;
                return 
            end
        end
    end
else
    dispv(1, '   No betas found, starting FIR1 computations...');
end

%%  Create target directory if not existing

if ~isfolder(target_dir_subj)
    dispv(1, ['   Creating target_dir_subj: ' target_dir_subj]);
    [s, ~] = mkdir(target_dir_subj); % s = succes (1)
    assert(s, 'Could not create target directory %s', target_dir_subj);
end

%% For FIR1: use duration 0 (todo: not a nice way to put it here)

duration = false; % set all durations to 0 if false

%% Get directories of the nifti files in the subject's folder

% Check for passed niftis (here, for normalized scans, can expand to other analysese)
if ~isfield(cfg, 'normalized_data') || isempty(cfg.normalized_data)
    
    % Check folder contents
    afx = cellstr(spm_select('FPList', source_dir_subj, '.*run-[0-9]{2}.nii'));
    assert(~isempty(afx), 'Could not recognize data in %s, please check.', source_dir_subj)
    dispv(1, '\n5. Found input data:'); disp(afx);
    
else
    assert(~isempty(cfg.normalized_data), ...
        'Could not recognize data in %s, please check.', cfg.normalized_data)
    afx = cfg.normalized_data;
end

%% Design setup

dispv(1, '6. Specifying design for sub-%02i...', subj)
dispv(2, '   Loading %s', tsv_file, 'display')
events = tdfread(tsv_file,'\t');

if ~duration
    dispv(2, '   Setting duration to 0')
    events.duration(:) = 0;
end

% Pre-allocate regressor information per run
design_files = cell(size(afx))'; % init
eventsize   = size(events.onset);
if eventsize(1) >= 2100 % rough estimation for now, but 7 runs should be <2100 betas
    runsize    = round(eventsize(1)/8);  
else
    runsize    = round(eventsize(1)/7); 
end

%% Get & order run information for SPM

for run_ind = 1:length(afx)
    
    curr_bname = spm_file(afx{run_ind}, 'basename');
    curr_fname = fullfile(target_dir_subj, [curr_bname,'_conditions-FIR1.mat']);

    %% Get information for current run (replaces old version)
    
    % init (see runsize pre-allocation)
    onsets  = cell(runsize, 1); durations   = cell(runsize, 1); 
    numbers = cell(runsize, 1); block_type  = cell(runsize, 1); 
    run     = cell(runsize, 1); difficulty  = cell(runsize, 1);
    period  = cell(runsize, 1); phase       = cell(runsize, 1); 
    names   = cell(runsize, 1); sort_var    = cell(runsize, 1);
    
    regr_ind = 0; % init regressor counter
    
    for event_ind = 1:size(events.onset, 1)
        
        % Get current data from tsv struct
        curr_event_run = events.run(event_ind);
        
        if curr_event_run ~= run_ind
            continue % current event from tsv file not for current run, skip it
        else
            % put current event as next event in SPM regressor info
            regr_ind = regr_ind + 1;
            
            % Pull events data from tsv struct into variables
            names{regr_ind} =   [events.block_type(event_ind, 1:3), ...
                                upper(events.difficulty(event_ind)), ...
                                events.period(event_ind, :), ...
                                sprintf('%.2i',events.phase(event_ind))];
            
            onsets{regr_ind}       = events.onset(event_ind);
            durations{regr_ind}    = events.duration(event_ind);
            numbers{regr_ind}      = events.number(event_ind);
            block_type{regr_ind}   = events.block_type(event_ind, 1:3);
            run{regr_ind}          = events.run(event_ind);
            difficulty{regr_ind}   = events.difficulty(event_ind);
            period{regr_ind}       = events.period(event_ind);
            phase{regr_ind}        = events.phase(event_ind);
            
            % Create map with values for each sorting conditions
            task_map        = {'Kon', 'Sum', 'Aeh';
                                10000, 20000, 30000};
            difficulty_map  = {'s', 'c';
                               1000, 2000}; % simple/complex
            period_map      = {'Cue  ', 'Trial', 'Rest ';
                               100 ,     200    , 300};
             
            % Get each of these values for the current event
            val_task       = task_map(2, strncmpi(task_map(1, :), ...
                names{regr_ind}, 3)); % matches Kon, Sum, Aeh   
            val_difficulty = difficulty_map(2, strncmpi(difficulty_map(1, :), ...
                names{regr_ind}(4), 1)); % matches S, C            
            val_period     = period_map(2, strncmpi(period_map(1, :), ...
                names{regr_ind}(5:7), 1)); % mathces cue, tri, res
            val_phase      = events.phase(regr_ind);
            
            % Sum all up (starts at 11100)
            sort_var{regr_ind}(1) = val_task{:} + val_difficulty{:} + val_period{:} + val_phase; 
        end
    end
    
    dispv(2, '  Sorting run info for sub-%02i...', subj)
    
    % Sort fields according to condition - diff - period - phase (phase
    % order is preserved)
    [~, sort_ind] = sort([sort_var{:}]);
    
    % sort 'names','onsets','durations','numbers'
    names       = names(sort_ind);
    onsets      = onsets(sort_ind);
    durations   = durations(sort_ind);
    numbers     = numbers(sort_ind);
    
    %% Write runwise
    
    dispv(2, ['   Writing designfile run ' num2str(run_ind) ' (' num2str(length(names)) ' regressors): ' curr_fname])
    save(curr_fname, 'names', 'onsets', 'durations', 'numbers', 'run');
    design_files{run_ind} = curr_fname;
    
end

dispv(1, '\nFirst-level design setup complete.\n')

%% Define directories for remote_setup (if selected)

if isfield(cfg, 'use_remote_setup_tmpdir')
    rs_cfg = {};
    rs_cfg.perm_basedir = target_dir_subj;
    % IMPORTANT!!! Make sure that the same base directory is not used by
    % different processes at the same time.
    % The following entry is optional. by default, it will be replaced with
    % /local/cfg.perm_basedir
    if ispc
        rs_cfg.tmp_basedir = fullfile(cfg.all_paths.localdir, ...
            strrep(rs_cfg.perm_basedir, ':\', '')); % e.g. c:\temp\...
    else
        rs_cfg.tmp_basedir = fullfile(cfg.all_paths.localdir, rs_cfg.perm_basedir);
    end
    % Start
    rs_cfg = remote_setup(rs_cfg); % Function made by Kai
end

%% Store input file information

clear json
json.func           = afx;
inputfiles_fname    = fullfile(target_dir_subj, 'input_files.json');
dispv(1, ['Writing ' inputfiles_fname])
spm_jsonwrite(inputfiles_fname, json, struct('indent', '\t'));

diary off % SPM output will probably be unreadable in diary

%% Design setup

clear matlabbatch

dispv(1, 'FIR regressor (beta) computation: Handing design to SPM...')
for run_ind = 1:length(afx)
    curr_fname = design_files{run_ind};
    matlabbatch{1}.spm.stats.fmri_spec.sess(run_ind).scans = cellstr(afx{run_ind}); % per run (loop)
    matlabbatch{1}.spm.stats.fmri_spec.sess(run_ind).multi = cellstr(curr_fname); % per run
    matlabbatch{1}.spm.stats.fmri_spec.sess(run_ind).hpf   = 10e6; % was: 128
end

matlabbatch{1}.spm.stats.fmri_spec.dir                 = cellstr(target_dir_subj);
matlabbatch{1}.spm.stats.fmri_spec.timing.units        = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT           = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t       = 33;  %warning('fmri_t and fmri_t0 still set to defaults (16 and 8); TODO: see if modification necessary due to slice time correction')
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0      = 16;  % reference slice?
matlabbatch{1}.spm.stats.fmri_spec.fact                = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length    = 2;
matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order     = 1;
matlabbatch{1}.spm.stats.fmri_spec.volt                = 1;
matlabbatch{1}.spm.stats.fmri_spec.global              = 'None';

% Set mask threshold to -Inf if you want to deactivate automatic masking
if isfield(cfg, 'saa') && isfield(cfg.saa, 'spmautomaticmasking') && ...
        cfg.saa.spmautomaticmasking == 1
    matlabbatch{1}.spm.stats.fmri_spec.mthresh         = 0; %.8;
elseif isfield(cfg, 'saa')
    matlabbatch{1}.spm.stats.fmri_spec.mthresh         = -Inf; %.8;
else 
    matlabbatch{1}.spm.stats.fmri_spec.mthresh         = 0.8;
end
matlabbatch{1}.spm.stats.fmri_spec.mask                = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi                 = 'None'; % was: AR(1)

%% Run batch
dispv(1, '   Saving design...')
spm_jobman('run', matlabbatch);
dispv(1, '   FIR design saved to %s --> SPM.mat', target_dir_subj);

%% Estimate

clear matlabbatch

matlabbatch{1}.spm.stats.fmri_est.spmmat            = cellstr(fullfile(target_dir_subj,'SPM.mat'));
matlabbatch{1}.spm.stats.fmri_est.write_residuals   = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1;

dispv(1, '   Estimate SPM mat...')
spm_jobman('run', matlabbatch);

dispv(1, '\nFirst level analysis completed.\n')

%% Copy mask files, sim_cfg and nii_nam.mat
% select files to copy

if strcmpi(cfg.data_type, 'simulated')
    cpfiles = {};
    cpfiles = [cpfiles; cellstr(spm_select('FPList', source_dir_subj, 'sub-.*mask.*.nii'))];
    cpfiles = [cpfiles; cellstr(spm_select('FPList', source_dir_subj, 'mask_names.mat'))];
    cpfiles = [cpfiles; cellstr(spm_select('FPList', source_dir_subj, 'nii_nam.mat'))];
    cpfiles = [cpfiles; cellstr(spm_select('FPList', source_dir_subj, 'sim_cfg.mat'))];
    
    % copy
    for cp_ind = 1:length(cpfiles)
        source = cpfiles{cp_ind};
        dispv(1, ['   Copying ' source ' --->\n   ' target_dir_subj]);
        copyfile(source, target_dir_subj);
    end
end

%% If tempo setup was selected, copy everything back
if isfield(cfg, 'use_remote_setup_tmpdir')
    rs_cfg = remote_teardown(rs_cfg);
end

end