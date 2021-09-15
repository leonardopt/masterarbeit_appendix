% function subj_cfg = generate_single_subject_data(cfg, sub_idx)
%
% Generates data per subject based on some value inputs for simulation
% and the general data structure of the empirical data.
%
% In:
% Subject-specific cfg (subj_cfg), data simulation parameters in
%  cfg.saa.data_gen
% Current subs_todo indexed subject (subj)
%
% Out:
% Full simulation cfg and simulated voxel data in nifti format with
% added noise and other parameters
%
% Needs:
% generate_time_series.m, assign_cue_voxel_weights.m,
% add_voxels_with_mixed_effects.m, make_SPM_compatible.m,
% generate_nifti_from_models.m, compute_sim_SNR.m,
% plot_generated_data.m

% History:
% Added SNR and plotting, changed some variable names -- Ingmar 03/09/20
% Made by Kai and Leonardo -- 2019
% Separated from generate_data -- Ingmar 31/03/20
% Separated from first level analysis functions -- Ingmar 08/04/20

function [subj_cfg, mask_names, nii_nam] = generate_single_subject_data(cfg, subj)
%% Unpack setting for model

n_mixedvoxels   = cfg.saa.data_gen.n_mixedvoxels;
weight          = cfg.saa.data_gen.weight_mixed_voxels;
models          = cfg.saa.data_gen.all_models;
modelnames      = cfg.saa.data_gen.model_names;
additionalnoise = cfg.saa.data_gen.additionalnoise;
prefweight      = cfg.saa.data_gen.pref_weight;

%% Create timeseries

[ts_data, regr_cue_names] = generate_time_series(cfg, subj);
dispv(1,'\nI. Generated voxel time series for subject %0i\n', subj');

%% Assign different weights to different voxels

% assign_voxel_weights also used a weightformat parameter which is no longer needed
weighted_data = assign_cue_voxel_weights(cfg, ts_data, regr_cue_names, models, ...
    modelnames, prefweight, additionalnoise);

dispv(1,'\nII. Assigned different weights to response voxels\n');

%% Add voxels with mixed effects

for d_ind = 1:length(weighted_data)
    weighted_data{d_ind} = add_voxels_with_mixed_effects(weighted_data{d_ind}, n_mixedvoxels(1), weight(d_ind));
end

dispv(1,'\nIII. Added mixed effects\n');

%% Make SPM compatible
% Make time-series compatible with SPM's weird way of calculating masks
% Do it after having added all voxels!

% Leos original version: Mask every layer
% if cfg.saa.spmautomaticmasking
%     for d_ind = 1:length(weighted_data)
%        error('Change here: change make_SPM_compatible to work with all layers (weighted_data) at once to get one common offset')
%        [weighted_data{d_ind}, make_SPM_compatible_info] = make_SPM_compatible(weighted_data{d_ind});
%     end
%    dispv(1, '\nIV. Design matrix made SPM compatible using an offset of %g\n', make_SPM_compatible_info.offset);
% end

[weighted_data, make_SPM_compatible_info] = make_SPM_compatible(weighted_data);

dispv(1, '\nIV. Design matrix made SPM compatible using an offset of %g\n', ...
    make_SPM_compatible_info.offset);

%% Create nifti images

% Generates nifti images for all subjects, on the basis of the cfg input models
[subj_cfg, mask_names, all_nii_dat, nii_nam] = generate_nifti_from_models(cfg, ...
    subj, weighted_data);

dispv(1, '\nV. Voxel data written in nifti images (or checked) for subject %0i.\n', subj);

%% Compute signal to noise ratio (SNR) for all runs/layers

subj_cfg = compute_sim_SNR(subj_cfg, weighted_data, make_SPM_compatible_info.offset);

dispv(1, '\nVI. Signal-to-noise ratio (SNR) computed with dB = snr() and sigmaS/sigmaN for sub-%02i.\n', subj);

%% Save

cfg             = subj_cfg; % Replace
preproc_sub     = fullfile(cfg.preproc_dir, sprintf('sub-%02i', subj));
cfg_file        = fullfile(preproc_sub, 'sim_cfg.mat');
dispv(1,'\n     Storing full sim_cfg to %s\n', cfg_file);
save(cfg_file, 'cfg');

%% Plot any required info

% If plotting is not needed
if isfield(cfg.saa.data_gen, 'plotting') && ~cfg.saa.data_gen.plotting
    dispv(1, '    Skipping subject plots since cfg.saa.data_gen.plotting == 0...\n');
    return % <--- Skip function
end

% Plots output provided here (new), but skip for larger loops >7 subs
if cfg.n_subs <= 7 || (cfg.n_subs > 7 && find(cfg.subs_todo == cfg.subj) < 7)
    plot_generated_data(subj_cfg, subj, ts_data, weighted_data, all_nii_dat, ...
        nii_nam, make_SPM_compatible_info.offset)
    dispv(1, '\nVII. Voxel data plotted');
    
% Log and skip    
else
    txtfile = fullfile(preproc_sub, 'subject_plots_only_in_first_7_sub_folders.txt');
    write_str = ['Individual plots were only generated for the first 7 ', ...
        'pipeline iterations because data generation parameters\n', ...
        'are equal across iterations and will thus vary little per "subject".\n', ...
        'Note: when running parfor, the skipped "subjects" might not', ...
        ' be in the original subject list order (after sub-07).\n\n'];
    fid = fopen(txtfile, 'a+');
    fprintf(fid, sprintf(write_str, subj));
    fclose(fid);
    % TODO: implement this for subject decoding plots, to avoid taking too
    % much space?
end

