% function [all_ACC, all_CONFMAT, first_data, masks_ind_todo, ...
% curr_ts_SNR] = afx_timeresolved_decoding_load_group_data(run_cfg, ...
% s_ind, sourcefolder, first_data, masks_ind_todo, mask_names)
%
% Loads all accuracy data for given subjects to do. Sorts data in structs
% per mask. Can be used with parfor for faster loading (since loading is
% slow in MATLAB). Will also try to load confusion matrices per subject.
%
% Kai, Jun 2020, adapted by Ingmar

function [all_ACC, all_CONFMAT, first_data, masks_ind_todo, curr_ts_SNR] = afx_timeresolved_decoding_load_group_data(run_cfg, ...
    s_ind, sourcefolder, first_data, masks_ind_todo, mask_names)
%% Verbosity, check loop size
% verbose seems to be lost/overwritten, dispv currently always displays;
% gives too much cmd output in grouplevel analyses. TODO: check why

% Attempt to force verbosity, not sure if this works or if somethings
% missing
if ~exist('verbose', 'var'), global verbose; verbose = 1; %#ok<TLEV>
elseif isfield(run_cfg, 'verbose'), verbose = run_cfg.verbose;
end

subs_todo = run_cfg.subs_todo;
if length(subs_todo) > 10 || ~run_cfg.group_template
    displvl = 2; % Reduce output to cmd (some strings are always set to 2)
else
    displvl = 1;
end

%% Get subjects, folders, cfg

curr_sub            = subs_todo(s_ind);
sourcefolder_sub    = fullfile(sourcefolder, sprintf('sub-%02i', curr_sub));
dispv(1, '    Loading data for subject %02i ...', curr_sub);
dispv(displvl, 'Loading from %s ...', sourcefolder_sub);

% Load cfg
curr_sbj_cfg_file   = fullfile(sourcefolder_sub, 'cfg.mat');
assert(isfile(curr_sbj_cfg_file), ...
    'Could not locate cfg.mat from %s. Check if decoding was completed correctly.', ...
    sourcefolder_sub);
curr_cfg = load(curr_sbj_cfg_file, 'cfg'); % <--- load
curr_cfg = curr_cfg.cfg;
dispv(2, '     Loaded cfg: %s', curr_sbj_cfg_file)

% Get masks of current subject; store masks to check later
full_mask_names = curr_cfg.files.mask;
if ~iscell(full_mask_names), full_mask_names = cell(full_mask_names); end % convert to cell

curr_mask_names = {};
for m_ind = 1:length(full_mask_names)
    [~, curr_mask_names{m_ind}, ~] = fileparts(full_mask_names{m_ind}); %#ok<*AGROW>
end

% Also remove subject specific part
curr_mask_names = strrep(curr_mask_names, sprintf('sub-%02i_', curr_sub), '');

if ~exist('mask_names', 'var')
    dispv(2, '     Retrieving mask names for all subjects from subject %02i...', curr_sub);
    mask_names = curr_mask_names; % as extracted above
else
    % Check masks
    assert(all(strcmp(mask_names, curr_mask_names)), ...
        'group_decoding:masknamesmatch', ...
        'Expected mask names do not match those found in source directory, sbj %i differs from sbj %i, check why.',  ...
        subs_todo(1), curr_sub);
end

% Get mask index if not present
if ~exist('masks_ind_todo', 'var')
    dispv(2, '     Retrieving masks_ind_todo from subject %02i...', curr_sub);
    masks_ind_todo = 1:length(curr_cfg.files.mask);
    dispv(displvl, '     Found %i masks.', length(masks_ind_todo));
end

%% Load all indicated masks

% Init
all_ACC         = [];
all_CONFMAT     = [];
curr_ts_SNR     = [];
curr_nii_SNR    = [];
curr_data       = {};
for m_ind = masks_ind_todo
    curr_maskname = mask_names{m_ind};
    curr_maskfile = ['res-ACC_' sprintf('sub-%02i', curr_sub) '_' curr_maskname '.mat'];
    fullmaskfile  = fullfile(sourcefolder_sub, curr_maskfile);
    
    if ~exist(fullmaskfile, 'file')
        dispv(2, '    Could not find %s.', fullmaskfile);
        curr_maskfile = ['res-ACC_' curr_maskname '.mat'];
        fullmaskfile  = fullfile(sourcefolder_sub, curr_maskfile);
        
        if ~exist(fullmaskfile, 'file')
            dispv(0, ['    Could not find filename or alternative without ', ...
                'sub-xx. Data might be missing. Skipping participant and logging error.']);
            all_ACC = []; all_CONFMAT = []; % to prevent output argument not assigned error
            
            try % Log in existing file
                if ~isfolder(fullfile(sourcefolder, 'group'))
                    mkdir(fullfile(sourcefolder, 'group'));
                end
                
                fname_errorlog = fullfile(sourcefolder, 'group', 'grouplevel_loading_failed.txt');
                write_str      = sprintf('%s: Subject %02i mask ACC data loading failed for group analysis', ...
                    datestr(now), curr_sub);
                fid = fopen(fname_errorlog, 'a+');
                fprintf(fid, ['\n' write_str]); % write on new line
                fclose(fid);
                dispv(1, '    Logged sub-%02i : %s', curr_sub, fname_errorlog);
            catch write_err % probably no writing permissions
                try
                    warning('Could not write errorlog to sourcedir, saving to ACC data')
                    disp(['    Error: ' write_err.message]);
                    all_ACC.error = write_str;
                catch
                    warning('Writing errorlog to sourcedir or ACC data failed.')
                end
            end
            
            return % < --- skip rest to try loading next sub
        end
        dispv(2, '    Alternative file name without sub-xx found.');
    end
    
    dispv(displvl, '  B. Loading masks...');
    dispv(2, '     Loading: %s', fullmaskfile);
    curr_data{m_ind} = load(fullmaskfile); % <--- load
    close all % Avoid that figures are loaded
    
    % for first sbj, init ACC
    if s_ind == 1, first_data = curr_data{m_ind}; end % keep as reference
    if s_ind == 1, first_cfg = curr_cfg; end
 
    %% Collect ACC images
    
    dispv(displvl, '  C. Building ACC images...');
    
    if curr_cfg.design.n_cond_per_step == 3 && ...
            curr_cfg.design.set2conditions_map(end, 1) == 1
        n_comps = 1;
    elseif curr_cfg.design.n_cond_per_step <= 3 % likely 2-step xclass or pairwise
        n_comps = length(curr_data{1}.ACC_mask_data.ACCs);
    end
    
    for a_ind = 1:n_comps
        dispv(displvl, '    Checking image %0i/%0i...', a_ind, n_comps)
        
        % Checks
        if ~isfield(run_cfg, 'disable_grouploading_checks') || ~run_cfg.disable_grouploading_checks
            if length(curr_data{m_ind}.ACC_mask_data.cfg) > 1 && s_ind > 1
                % Basic check that dimensions are the same: not super relevant
                % though... all data should be 45 timepoints and will have the
                % same timepoint labels coming from make_design_timeresolved
                assert(isequal(curr_cfg.design.conditions2set_map_axislabel{a_ind}.traintime_yaxis_str, ... % yaxis str
                    first_cfg.design.conditions2set_map_axislabel{a_ind}.traintime_yaxis_str), ...
                    'Labels cfg.design.conditions2set_map_axislabel{a_ind}.traintime_yaxis_str disagree between sbj %i and %i', ...
                    subs_todo(1), subs_todo(s_ind));
                assert(isequal(curr_cfg.design.conditions2set_map_axislabel{a_ind}.testtime_xaxis_str, ... % xaxis str
                    first_cfg.design.conditions2set_map_axislabel{a_ind}.testtime_xaxis_str), ...
                    'Labels cfg.design.conditions2set_map_axislabel{a_ind}.testtime_xaxis_str disagree between sbj %i and %i', ...
                    subs_todo(1), subs_todo(s_ind));
            end
            
            assert(isequal(curr_cfg.design.conditions_per_step, ...
                first_cfg.design.conditions_per_step), ...
                'Difference cfg.design.conditions_per_step between sbj %i and %i', ...
                subs_todo(1), subs_todo(s_ind));
            assert(isequal(curr_cfg.design.label_val2str_map, ...
                first_cfg.design.label_val2str_map), ...
                'Difference cfg.design.conditions_per_step between sbj %i and %i', ...
                subs_todo(1), subs_todo(s_ind));
        end % checks
        
        % Get current ACC result
        if isstruct(curr_data{m_ind}.ACC_mask_data)
            curr_ACC = curr_data{m_ind}.ACC_mask_data.ACCs{a_ind};
        else
            curr_ACC = curr_data{m_ind}.ACC_mask_data{a_ind}.ACCs{a_ind};
        end
        
        % Add data
        all_ACC{m_ind}.ACCs{a_ind}.ACC(:, :, s_ind) = curr_ACC.ACC;
        
        % Set/check title
        if ~isfield(all_ACC{m_ind}.ACCs{a_ind}, 'title')
            all_ACC{m_ind}.ACCs{a_ind}.title = curr_ACC.title;
        end
        assert(isequal(all_ACC{m_ind}.ACCs{a_ind}.title, curr_ACC.title), ...
            'all_ACC{m_ind}.ACCs{a_ind}.title =/= curr_ACC.title')
        if ~isfield(all_ACC{m_ind}.ACCs{a_ind}, 'axislabel')
            all_ACC{m_ind}.ACCs{a_ind}.axislabel = curr_ACC.axislabel;
        end
        assert(isequal(all_ACC{m_ind}.ACCs{a_ind}.axislabel, curr_ACC.axislabel), ...
            'all_ACC{m_ind}.ACCs{a_ind}.axislabel =/= curr_ACC.axislabel')
    end
    
    if isfield(all_ACC{m_ind}, 'mask_name')
        assert(isequal(all_ACC{m_ind}.mask_name, curr_maskname), ...
            'Unexpected: masknames disagree (%s vs %s), please check why', ...
            all_ACC{m_ind}.mask_name, curr_maskname);
    else
        all_ACC{m_ind}.mask_name = curr_maskname;
    end
    
    dispv(displvl, '    Loaded data passed basic checks...');
    all_ACC{m_ind}.full_datafile{s_ind} = fullmaskfile; % Place in cell structs
end

%% Confusion matrix loading,...

% We don't need extensive checks because the ACCs were already checked.
% Means script will preload one sub for checking the parfor
if isfield(run_cfg, 'group_template') && run_cfg.group_template
    all_CONFMAT = []; curr_ts_SNR = [];
    return % can skip
end

if (any(strcmpi(run_cfg.results.output, 'confusion_matrix')) || ...
        any(strcmpi(run_cfg.results.output, 'confusion_matrix_plus_undecided'))) && ...
        ((isfield(run_cfg.design, 'make_pairwise') && run_cfg.design.make_pairwise == 0) || ...
        ~isfield(run_cfg.design, 'make_pairwise')) 
    for m_ind = masks_ind_todo
        curr_maskname = mask_names{m_ind};
        curr_maskfile = ['res-CONFMAT_' sprintf('sub-%02i', curr_sub) '_' curr_maskname '.mat'];
        confmat_maskfile  = fullfile(sourcefolder_sub, curr_maskfile);
        
        if ~exist(confmat_maskfile, 'file')
            dispv(2, '  Could not find ...%s.', confmat_maskfile);
            curr_maskfile = ['res-CONFMAT_' curr_maskname '.mat'];
            confmat_maskfile  = fullfile(sourcefolder_sub, curr_maskfile);
            
            if ~exist(confmat_maskfile, 'file')
                dispv(0, [ '  Could not find CONFMAT file or alternative without ', ...
                    'sub-xx in name. Data might be missing. Skipping participant and logging error']);
                all_CONFMAT = [];
                
                try % Log in existing file
                    fname_errorlog = fullfile(sourcefolder, 'group', 'grouplevel_loading_failed.txt');
                    write_str      = sprintf('%s: Subject %02i mask CONFMAT data loading failed for group analysis', ...
                        datestr(now), curr_sub);
                    fid = fopen(fname_errorlog, 'a+');
                    fprintf(fid, ['\n' write_str]); % write on new line
                    fclose(fid);
                    dispv(1, '   Logged sub-%02i : %s', curr_sub, fname_errorlog);
                catch write_err % probably no writing permissions
                    disp(['    Error: ' write_err.message]);
                    try
                        all_CONFMAT.error = write_str; % might not be accessible to parfor
                        warning('Could not write errorlog to sourcedir, saving to ACC data')
                    catch
                        warning('Writing errorlog to sourcedir or ACC data failed.')
                    end
                end
                return % <--- skip rest to try loading next sub
            else
                dispv(2, '  Alternative file name without sub-xx found.');
            end
        end
        
        dispv(displvl, '    Loading confusion matrices for mask %i (%s)...', ...
            m_ind, curr_maskname);
        dispv(2, '   Loading: %s', confmat_maskfile);
        all_CONFMAT{m_ind} = load(confmat_maskfile);
        close all % Avoid that figures are loaded
        
    end % masks
else
    all_CONFMAT = []; 
end

%% SNR loading

% Check if any SNR data present
if strcmpi(curr_cfg.data_type, 'simulated') && isfield(curr_cfg.saa, ...
        'snr_db_total') && ~isempty(curr_cfg.saa.snr_db_total)
    try
        dispv(displvl, '    SNR values for simulated data available in curr_sbj_cfg, retrieving...');
        
        % Get total mean SNR in dB
        snr_db_total = curr_cfg.saa.snr_db_total;
                
        % Get total CNR
        CNR_total = curr_cfg.saa.CNR_total;
                
        % Format explicitly for parfor (struct did not work)
        % DO NOT CHANGE ORDER OF VARIABLES
        curr_ts_SNR = {snr_db_total, CNR_total};
        
    catch e1
        dispv(displvl, '    SNR value loading failed, see error:');
        disp(e1.message);
        curr_ts_SNR = []; % return empty
    end
end

% try this for "output argument" error?
if ~exist('curr_ts_SNR', 'var') || isempty(curr_ts_SNR)
    curr_ts_SNR = [];
end
if ~exist('curr_nii_SNR', 'var') || isempty(curr_nii_SNR)
    curr_nii_SNR = [];
end