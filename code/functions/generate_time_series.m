% function [ts_data, regr_cue_names] = generate_time_series(cfg, sub_idx)
%
% Generates time series for simulated voxels. Adds noise specified
% in cfg.saa.data_gen parameters
%
% IN:
%  cfg with the following subfields: saa.settings.conditions/subs_todo,
%   and saa.data_gen
%  subj : current subject in processing loop
%
% OUT:
%  data (formerly dat): structure which comprises every time series generated from
%   make_time_series_model_from_cue.m, which convolves raw timepoints with a
%   canonical HRF with added noise pre- and post- HRF.
%  regr_names: regressor names of the generated voxel time series
%
% Explanation and edits by Ingmar 04/2020
% Made by Kai/Leonardo, ~2019

function [all_ts, regr_cue_names] = generate_time_series(cfg, subj)
%% Get some variables from the cfg struct

tasktypes           = cfg.saa.tasks;
models              = cfg.saa.data_gen.all_models;
modelnames          = cfg.saa.data_gen.model_names;
hrf_conv            = cfg.saa.data_gen.hrf_conv;
noise_before_conv   = cfg.saa.data_gen.noise_before_conv;
noise_after_conv    = cfg.saa.data_gen.noise_after_conv;
layer_number        = cfg.saa.data_gen.layer_number;

%% Check that no model has nans (and warn if infs)

dispv(1, '\n    Checking model properties...');
for model_ind = 1:length(models)
    % Check that no model has nans and warn if infs
    assert(all(~isnan(models{model_ind})), ...
        'Found nans in model{%i}, please check why', model_ind);
    assert(any(~isinf(models{model_ind})), ...
        'Found infs in model{%i}, no idea why. Please check if it works as expected', model_ind);
    assert(~isempty(models{model_ind}(:)), ...
        'Model %g is empty, check pipeline', models{model_ind});
end

%% Retrieve SPM.mat (do not save to folder, it takes too much space)

spm_mat_dir     = cfg.level1_fmrispmmat_dir;
sub_str         = sprintf('sub-%02i', subj); % Set target directory for current subject
SPM_mat_source  = fullfile(spm_mat_dir, sub_str, 'func', 'SPM.mat'); % THIS IS THE SPM.mat FROM THE ORIGNAL fMRI FIR ANALYSIS

dispv(1, '    Loading %s ...', SPM_mat_source);

loaded_SPM  = load(SPM_mat_source);
SPM         = loaded_SPM.SPM;

dispv(1, '    Loading SPM.mat done\n');

%% Make time series for current participant

timeseries = cell(1, length(models)); % init

dispv(2, '    Creating timeseries from SPM.mat...')

% Loop for batching multiple noise parameters per layer in one nifti
for layer_ind = 1:length(layer_number) % length of parameter inputs
    
    for m_ind = 1:length(models)
        % ------------------------------------------------------------------- %
        timeseries{1, m_ind} = make_time_series_model_from_cue(cfg, SPM, tasktypes, models(m_ind), modelnames{m_ind}, hrf_conv, noise_before_conv(layer_ind), noise_after_conv(layer_ind));
        % ------------------------------------------------------------------- %
        regr_cue_names{m_ind} = timeseries{m_ind}(1).regr_cue_names; %#ok<*AGROW>
        regr_cue_names{m_ind} = strcat(regr_cue_names{m_ind}, '-', modelnames{m_ind});
    end
        
    %% Change data format
    
    for ts_ind = 1:length(timeseries{1})
        
        ts_data(ts_ind).X = []; % init
        for data_ind = 1:length(timeseries)
            
            n_regs = size(ts_data(ts_ind).X, 2); % preallocate
            curr_X = timeseries{data_ind}(ts_ind).X; % init
            
            for reg_ind = 1:size(curr_X, 2)
                
                % Concatenate model data
                ts_data(ts_ind).X(:, n_regs + reg_ind) = [curr_X(:, reg_ind)];
                
                % Concat cue regressor names
                ts_data(ts_ind).regr_cue_names(:, n_regs + reg_ind) = ...
                    timeseries{data_ind}(ts_ind).regr_cue_names(reg_ind);
            end
        end
    end
    
    % Check data validity
    dispv(2, '    Validating data (not empty, zero or NaN)...')
    assert(~isempty(ts_data),  'Data input is empty! Check pipeline');
    assert(~all(isnan(ts_data(ts_ind).X), 'all'), 'Data contains NaNs! Check pipeline input.');
    assert(~all(ts_data(ts_ind).X == 0, 'all'), 'Data is only zero! Check pipeline input.');
    
    all_ts{layer_ind} = ts_data; % pack in cell array
    
end % all ts

