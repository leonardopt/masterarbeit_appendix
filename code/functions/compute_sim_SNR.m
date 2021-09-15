% function cfg = compute_sim_SNR(cfg, weighted_data, offset)
%
% Computes mean SNR per run of the given generated data. Removes SPM global
% offset to zero-mean data.
%
% IN: weighted_data (timeseries data), nii_dat (nifti data), cfg.saa.data_gen
% OUT: cfg.saa.SNR, cfg.saa.SNR_CNR_descr
%
% Ingmar, 03-09-20

% Update by Ingmar to calculate zero-centered SNR and CNR for weighted 
%  data, excluding the cue, 08-01-21

function cfg = compute_sim_SNR(cfg, weighted_data, offset)
%% Options

do.noisemean = 0; % reduce noise variance with mean over all runs

%% Get parameters

layer_number    = cfg.saa.data_gen.layer_number;
layer_descr     = cfg.saa.data_gen.layer_descr;

dispv(1, ['    Computing signal reference from first run, first layer model', ...
    ' (assuming first model is noise-free!)...']);

% get minimum
for l_ind = 1:length(layer_number) % ? layers
    for r_ind = 1:length(weighted_data{1}) % 8 runs
        min_runlength(l_ind, r_ind) = size(weighted_data{l_ind}(r_ind).X, 1); 
    end
end
min_runlength = min(min(min_runlength));

%% Get timeseries SNR for response voxels

% pre-allocate
snr_db = zeros(1, size(layer_number, 2)); CNR = zeros(size(snr_db));
ignore_cue = 4; last_response_voxel = 9;

for l_ind = 1:length(layer_number) % ? layers
    
    % Mean over runs
    for r_ind = 1:length(weighted_data{1}) % 8 runs
        vox_modeldata(:, :, r_ind)  = weighted_data{1}(r_ind).X(1:min_runlength, 1:9) - round(offset); %#ok<*AGROW>
        vox_noise_data(:, :, r_ind) = weighted_data{l_ind}(r_ind).X(1:min_runlength, 1:9) - round(offset); % all runs
        vox_noise_est(:, :, r_ind)  = vox_noise_data(:, :, r_ind) - vox_modeldata(:, :, r_ind);
    end
    
    % Compute runwise and voxelwise
    if do.noisemean
        vox_mean_noise_est = mean(vox_noise_est, 3);
        for v_ind = ignore_cue:last_response_voxel % exclude cue
            for r_ind = 1:length(weighted_data{1}) % 8 runs
                vox_SNR(v_ind) = snr(vox_modeldata(:, v_ind, r_ind), vox_mean_noise_est(:, v_ind));
                vox_CNR(v_ind) = std(vox_modeldata(:, v_ind, r_ind) / std(vox_mean_noise_est(:, v_ind)));
            end
        end
    else
        for v_ind = ignore_cue:last_response_voxel % exclude cue
            for r_ind = 1:length(weighted_data{1}) % 8 runs
                vox_SNR(v_ind) = snr(vox_modeldata(:, v_ind, r_ind), vox_noise_est(:, v_ind, r_ind));
                vox_CNR(v_ind) = std(vox_modeldata(:, v_ind, r_ind) / std(vox_noise_est(:, v_ind, r_ind)));
            end
        end
    end
    
    % average over voxels
    snr_db(l_ind)  = mean(vox_SNR(ignore_cue:end), 'omitnan');
    CNR(l_ind)     = mean(vox_CNR(ignore_cue:end), 'omitnan');
end % layers

% Get total mean SNR in dB and CNR
cfg.saa.snr_db_total = snr_db;
cfg.saa.CNR_total    = CNR;

%% Display

dispv(1, '\n    SNR in dB for model (layer) 1 to %i:', length(layer_number));
disp(layer_descr); disp([cfg.saa.snr_db_total]);
dispv(1, '\n    CNR in sigmaS/sigmaN for model (layer) 1 to %i:', length(layer_number));
disp(layer_descr); disp([cfg.saa.CNR_total]);

pause(1);
