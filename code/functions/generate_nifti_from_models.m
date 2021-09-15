% function [cfg, mask_names, all_nii_dat, nii_nam] = generate_nifti_from_models(cfg, subj, all_dats)
%
% Function that generates nifti images for one subject, given the cell of
% timeseries all_dats. Computation will be skipped if folder already exists.
%
% IN: cfg, subject number, all_dats (timeseries data from model generation)
%
% OUT: cfg, mask_names, all_nii_dat (nifti data from transformed timeseries
%      data, nii_nam (names of each corresponding voxel location in 2d)
%
% requires: SPM, "create_nifti.m", "get_time_series.m"
%
% Edited by Ingmar, 04/05/2020
% Created by Kai/Leonardo, ~2019

function [cfg, mask_names, all_nii_dat, nii_nam] = generate_nifti_from_models(cfg, subj, all_dats)
%% Get some settings from the cfg struct, add folders if not there

sim_name        = cfg.sim_name;
preproc_dir     = cfg.preproc_dir;
preproc_sub     = fullfile(preproc_dir, sprintf('sub-%02i', subj));

%% Create directory for current condition or skip computation

% create data if flagged to overwrite
if isfield(cfg, 'flags') && isfield(cfg.flags, 'data_gen_overwrite_niftis') && cfg.flags.data_gen_overwrite_niftis
    dispv(1,'    Overwriting nii files, because cfg.flags.data_gen_overwrite_niftis = true')
    
    % try to load data
else
    try
        %% Loading
        disp(['    Found existing .nii, attempting to skip .nii file generation ', ...
            '(and discard currently generated data). ',...
            'To avoid this, set cfg.flags.data_gen_overwrite_niftis = true'])
        
        mask_names_fname = fullfile(preproc_sub, 'mask_names.mat');
        if exist(mask_names_fname, 'file')
            disp(['    Loading ' mask_names_fname])
            load(mask_names_fname, 'mask_names', 'nii_model_descr');
        else
            error('Could not load mask names, abort loading')
        end
        
        %% Checks
        
        nii_nam_fname = fullfile(preproc_sub, 'nii_nam.mat');
        if exist(nii_nam_fname, 'file')
            disp(['    Loading ' nii_nam_fname])
            load_nn = load(nii_nam_fname, 'nii_nam', 'nii_model_descr');
            assert(isequal(nii_model_descr, load_nn.nii_model_descr), ...
                'nii_model_descr does not agree between nii_nam.mat and mask_names.mat. Recomputing')
            nii_nam = load_nn.nii_nam;
        else
            error('Could not load mask names, abort loading')
        end
        
        stored_cfg = fullfile(preproc_sub, 'sim_cfg.mat');
        if exist(stored_cfg, 'file')
            disp(['    Loading ' stored_cfg])
            loaded_cfg = load(stored_cfg, 'cfg');
        else
            error('Could not load sim cfg, abort loading')
        end
        
        dispv(1, '\n    Checking data...')
        
        assert(isequal(cfg.data_type, loaded_cfg.cfg.data_type), ...
            'isequal(cfg.data_type, loaded_cfg.cfg.data_type) not met. Files might be corrupt, recomputing.')
        assert(isequal(cfg.sim_name, loaded_cfg.cfg.sim_name), ...
            'isequal(cfg.sim_name, loaded_cfg.cfg.sim_name) not met. Files might be corrupt, recomputing.')
        assert(isequal(nii_model_descr, loaded_cfg.cfg.nii_model_descr), ...
            'nii_model_descr did not agree between mask_names.mat and loaded cfg, recomputing');
        assert(isequal(mask_names, loaded_cfg.cfg.mask_names), ...
            'mask_names did not agree between mask_names.mat and loaded cfg, recomputing');
        if exist('all_dats', 'var') % Input can be missing
            datacheck = all_dats{1}.X;
            assert(~isempty(datacheck),  'all_dats{1}.X input is empty! Check pipeline');
        end
        
        dispv(1, '    Data validated.')
        
        try % Not crucial, some info might be missing
            dispv(1, '    Augmenting cfg with all information added during creation...')
            cfg.mask_names = mask_names; % checked above
            cfg.nii_model_descr = loaded_cfg.cfg.nii_model_descr;
            loaded_cfg.cfg.general_info.n_voxels_before_spm_masking = size(all_dats{1}(1).X, 2);
            cfg.general_info.n_voxels_before_spm_masking = loaded_cfg.cfg.general_info.n_voxels_before_spm_masking;
            
            % add info that data was loaded from cache
            cfg.generated_data_loaded_time_log = loaded_cfg.cfg.time_log;
        catch
            dispv(1, '    Not all info was succesfully added to cfg. Check later.');
        end
        
        all_nii_dat = {}; dispv(1, '    nii_dat run data not available');
        
        return % <---- Loading finished, skipping rest of function
        
    catch errloadnii
        warning('Loading nii files aborted due to error, generating new ni.')
        disp(['    Error: ' errloadnii.message]);
        clear mask_names
    end
end

%% Define directories for remote_setup (if selected)

if isfield(cfg, 'use_remote_setup_tmpdir')
    rs_cfg = {}; % init
    rs_cfg.perm_basedir = preproc_sub;
    
    % IMPORTANT! Make sure that the same base directory is not used by
    % different processes at the same time. The following entry is optional.
    % By default, it will be replaced with /local/cfg.perm_basedir
    if ispc
        rs_cfg.tmp_basedir = fullfile(cfg.use_remote_setup_tmpdir, strrep(rs_cfg.perm_basedir, ':\', '')); % remove second :\ if present
    else
        rs_cfg.tmp_basedir = fullfile(cfg.use_remote_setup_tmpdir, rs_cfg.perm_basedir);
    end
    
    preproc_sub = rs_cfg.tmp_basedir; % will be created, if it does not exist, and be deleted at the end
    rs_cfg      = remote_setup(rs_cfg); % Start
end

%% Create nifti file for each run

dispv(1,'\n    A. Create nifti files for each run...\n');
alldat_merged   = [all_dats{:}];
n_runs          = length(alldat_merged)/length(all_dats);
all_nii_dat     = cell(size(n_runs));

for run_number = 1:n_runs
    dispv(2,'    Run %i processing...', run_number);
    
    curr_nifti      = sprintf('sub-%02i_%s_run-%02i', subj, sim_name, run_number);
    fname           = fullfile(preproc_sub, curr_nifti);
    dats_to_write   = alldat_merged(run_number:n_runs:length(alldat_merged));
    
    % Creates nifti image
    % ------------------------------------------------------------------- %
    [nii_dat, nii_nam]  = create_nifti_multiple_models(dats_to_write, fname);
    % ------------------------------------------------------------------- %
    dispv(2,'    Run %i done.', run_number);
    
    % Save for computing SNR
    all_nii_dat{run_number} = nii_dat;
end
dispv(1,'    Nifti images for subject %02i - created.\n', subj);

%% Create mask

dispv(1, '\n    B. Generating masks for current subject...\n');
% Separated from current script for easy access - 11/06
% ------------------------------------------------------------------- %
[cfg, mask_names] = create_masks(cfg, nii_nam, nii_dat, subj);
% ------------------------------------------------------------------- %
dispv(1, '    Generating all masks for current subject... done\n');

%% Add early info on full timepoints to cfg

dispv(1, '\n    C. Adding some info to cfg and json...\n')

% all timepoints minus constants minus 5 possible "missing" betas (276-271) per run
cfg.saa.data_gen.n_runs              = n_runs;
cfg.saa.data_gen.t_runlength         = size(nii_dat, 4);
cfg.saa.data_gen.expected_betacount  = (n_runs * size(nii_dat, 4)) - (1*n_runs) - (5*n_runs);

%% save: README.txt that describes all files/types

% add to info
info.sim_cfg = 'cfg with parameters used in generate_nifti_from_models.m';
info_fname = fullfile(preproc_sub, 'info.txt');
disp(['    Saving info to ' info_fname])
info_file = fopen(info_fname, 'w');
fwrite(info_file, jsonencode(info));
fclose(info_file);

%% If remote (temporary) setup was selected, copy everything back

if isfield(cfg, 'use_remote_setup_tmpdir')
    rs_cfg = remote_teardown(rs_cfg);
end

%% Done

dispv(1, '    Finished nifti generation from models for subject %02i.', subj);

end