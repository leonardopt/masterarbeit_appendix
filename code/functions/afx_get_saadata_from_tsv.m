% afx_get_saadata_from_tsv()

function afx_get_saadata_from_tsv(curr_sbj, TR, preproc_sub, betas_per_sess)

% set path to tdt-saa
global all_paths
global_all_paths
set_all_paths
% init saa_cfg
saa_cfg = saa_defaults;

%% Set parameters
% set result directory
% preproc_sub = fullfile(all_paths.data_base_dir, sprintf('derivatives_test/preproc-saatsv/sub-%02i/func/', curr_sbj))

%% load design data
tsv_filename = fullfile(all_paths.data_base_dir, sprintf('sub-%02i/func/sub-%02i_task-afx_events.tsv', curr_sbj, curr_sbj));
tsv_data = saa_load_data(tsv_filename); 

%% load & apply testsuitedef
suitedefdir = fileparts(which('afx_get_saadata_from_tsv.m')); % dir of this function
suitedeffile = fullfile(suitedefdir, 'afx_get_saadata_suitedef.tsv');
saa_cfg.testsuitedef = saa_load_testsuitedefinition(suitedeffile);

% apply testsuite and check entries
[tsv_data, saa_cfg] = saa_apply_testsuitedef(tsv_data, saa_cfg); % set trial and session var
saa_cfg = saa_checkvars(saa_cfg, tsv_data, {'session_var', 'trial_var'});

%% split tsv_data for trial and rest part
tsv_trial = tsv_data(~strcmp(tsv_data.onset, 'n/a'), :);
% tsv_rest = tsv_data(strcmp(tsv_data.onset, 'n/a'), :);
saa_cfg_trial = saa_cfg;
% saa_cfg_rest = saa_cfg;

%% augment data
% % augment data: full tsv (trial and rest mixed)
% [tsv_data, saa_cfg] = saa_dummyvar(tsv_data, saa_cfg);
% tsv_data = saa_add_prev(tsv_data, saa_cfg);
% % tsv_data = saa_add_synthetic_data(tsv_data, saa_cfg);
% % get testcases, will be used to create masks
% saa_cfg = saa_get_default_testcase_definitions(tsv_data, saa_cfg);

% augment trial data
[tsv_trial, saa_cfg_trial] = saa_dummyvar(tsv_trial, saa_cfg_trial);
tsv_trial = saa_add_prev(tsv_trial, saa_cfg_trial);
% tsv_trial = saa_add_synthetic_data(tsv_trial, saa_cfg_trial);
% get testcases, will be used to create masks
saa_cfg_trial = saa_get_default_testcase_definitions(tsv_trial, saa_cfg_trial);

% % augment rest data
% [tsv_rest, saa_cfg_rest] = saa_dummyvar(tsv_rest, saa_cfg_rest);
% tsv_rest = saa_add_prev(tsv_rest, saa_cfg_rest);
% % tsv_rest = saa_add_synthetic_data(tsv_rest, saa_cfg_rest);
% % get testcases, will be used to create masks
% saa_cfg_rest = saa_get_default_testcase_definitions(tsv_rest, saa_cfg_rest);

%% Create timecourses
warning('Todo: change vector with actual number of session for this subject (as vector betas_per_sess(sess_ind)=nimgs)')
% betas_per_sess = 276;
[timeresolved_data, labelnames, ext_labelnames, masks, saa_cfg] = afx_saatable2timecourses(saa_cfg_trial, tsv_trial, TR, betas_per_sess);
% timeresolved_data: n_betas x labelnames, labelnames: labels cellstr

%% %     % visualise (if you like)
% %% plot (maybe with color limits)
% curr_sess = 1
% figure('name', sprintf('sess %i (color limits)', curr_sess), 'position', [1           1        1600         800])
% pX = timeresolved_data{curr_sess};
% pL = labelnames;
% imagesc(pX);
% title(sprintf('sess %i', curr_sess));
% ylabel(sprintf('beta image (%i)', size(pX, 1)));
% set(gca, 'XTick', 1:size(pX, 2), 'XTickLabel', pL, 'XTickLabelRotation', 90);
% set(findall(gcf,'-property','FontSize'),'FontSize', 6)
% set(gca, 'clim', quantile(pX(:), [.3 .97]))
% colorbar
% drawnow

%% This should then give us our time series model + the masks: 
% continue with the functions that we already have

% run and create the model, store whatever is needed in a cache mat file to
% see what the exact timeseries size is (per run i guess). then try to
% integrate this into the pipeline.

% load('all_dats_cache.mat'); % example data for tim series to write, see
% start

% put into all_dats
all_dats = [];
for sess_ind = 1:length(timeresolved_data)
    all_dats(sess_ind).X = timeresolved_data{sess_ind};
    all_dats(sess_ind).voxelnames = ext_labelnames;
end

sim_name = 'saa-core';
if ~exist(preproc_sub, 'dir'), mkdir(preproc_sub); end

for sess_ind = 1:length(timeresolved_data)
    dispv(2,'    Run %i processing...', sess_ind);
    curr_nifti      = sprintf('sub-%02i_%s_run-%02i', curr_sbj, sim_name, sess_ind);
    fname           = fullfile(preproc_sub, curr_nifti);
    dats_to_write   = all_dats(sess_ind); % changed from multiple models, only 1 session
    
    % Creates nifti image
    % ------------------------------------------------------------------- %
    [nii_dat, nii_nam]  = create_nifti_multiple_models(dats_to_write, fname, false);
    % ------------------------------------------------------------------- %
    dispv(2,'    Run %i done.', sess_ind);
end

%% Create masks, write them as .nii and add maskname.mat
maskdim = size(nii_dat(:, :, :, 1));
vox_names = nii_nam(:, :, :, 1);
vox_names(cellfun('isempty', vox_names)) = {''};

dispv(2,'    Writing masks...');
mask_names = {};
for m_ind = 1:length(masks)
    curr_mask = masks(m_ind);
    [~, mask_idx] = intersect(vox_names, curr_mask.ext_labelnames);
    mask_vol = false(maskdim);
    mask_vol(mask_idx) = true;
    assert(isequal(unique(curr_mask.ext_labelnames), unique(vox_names(mask_vol)')), 'Unexpected error: masking the voxel names with mask_vol does not produce the voxel name(s) that were used to generate the mask. That should not happen');
    description = ['labels:', sprintf('%s;', curr_mask.labelnames{:}), ...
        '|ext_labels:' sprintf('%s;', curr_mask.ext_labelnames{:}), ...
        '|description:' curr_mask.descrip, '. ' ...
        '|expectation:' curr_mask.expect,  '. ' ...
        '|' curr_mask.info];

    curr_nifti = sprintf('sub-%02i_%s_mask-%i_name-%s.nii', curr_sbj, sim_name, m_ind, matlab.lang.makeValidName(curr_mask.maskname));
    fname      = fullfile(preproc_sub, curr_nifti);
    mask_names{end+1} = curr_nifti; %#ok<AGROW>
    %% Create data
    % Create image
    ni = nifti;

    % header info
    ni.dat       = file_array(fname, size(mask_vol), [0 spm_platform('bigend')], 0, 1, 0); % SPM
    ni.mat       = eye(4);
    ni.descrip   = description;
    ni.dat.dtype = 'FLOAT32';

    dispv(2, '    Writing: %s', fname);
    % prepare
    create(ni); % SPM
    % write image
    ni.dat(:,:,:) = mask_vol(:,:,:);
    if m_ind == 1
        full_mask_vol = mask_vol;
    else
        full_mask_vol = mask_vol | full_mask_vol; % add new mask to all masks
    end
end

% write mask.nii (full mask)
description_mask    = 'mask for all voxels of any simulation';
fname_mask          = fullfile(preproc_sub, 'mask.nii');
ni = nifti;
ni.dat       = file_array(fname_mask, size(full_mask_vol), [0 spm_platform('bigend')], 0, 1, 0); % SPM
ni.mat       = eye(4);
ni.descrip   = description_mask;
ni.dat.dtype = 'FLOAT32';
create(ni); % SPM
ni.dat(:,:,:) = full_mask_vol(:,:,:); % write image

% write masknames.mat & metadata (names of voxels)
masknames_fname = fullfile(preproc_sub, 'mask_names.mat');
nii_model_descr = vox_names;
disp(['    Saving names of masks for different models to ' masknames_fname])
info  = [];
info.mask_names         = 'Name of maskfiles for different simulations';
info.nii_model_descr    = 'Names of simulation in each voxel (across time, first vol in nii_nam)';
save(masknames_fname, 'mask_names', 'nii_model_descr', 'info');

% save description of each voxel
info.nii_nam = 'Names of simulation in each voxel (nii_model_descr)';
niinam_fname = fullfile(preproc_sub, 'nii_nam.mat');
disp(['    Saving  description of each voxel (nii_model_descr) to ' niinam_fname])
save(niinam_fname, 'nii_model_descr');

dispv(2,'    Writing masks done. afx_get_saadata_from_tsv done');
dispv(2,'afx_get_saadata_from_tsv done');
