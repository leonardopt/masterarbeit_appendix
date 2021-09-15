% function all_data = assign_cue_voxel_weights(cfg, ts_data, regnames, models, modelnames, prefweight, additionalnoise)
%
% Function to assign weight to voxels based on cue regressors returned from
% make_time_series_model_from_cue.m --> generate_time_series.m with some
% data structure reworking (based on SPM formatting which is quite cryptic)
%
% IN: data, regnames, models, modelnames
% OPTIONAL: prefweight, additionalnoise,
%  cfg.saa.data_gen.partialresponsevoxels == 1/0 (let voxels selective to a 
%   task respond to the other task signals or not)
% OUT: all_data
%
% TODO: expand at some point to include difficulty levels.. here?
%
% Kai/Leonardo ~2019, tweaked by Ingmar, June 2020

% Reverted back to February git version because partialresponsevoxels is a
% setting that works better for decoding. It makes the voxel selective to a
% certain condition not respond to other conditions, perhaps by (?)
% improving the discriminibality of classes during decoding. This in turn
% allows using higher noise levels etc. - Ingmar 21-07-21

function all_data = assign_cue_voxel_weights(cfg, all_ts, regr_cue_names, ...
    models, modelnames, prefweight, additionalnoise)
%% Development

% Default = 1;
if ~isfield(cfg.saa.data_gen, 'partialresponsevoxels') || isempty(cfg.saa.data_gen.partialresponsevoxels)
    cfg.saa.data_gen.partialresponsevoxels = 1; % set
end

%% Set

task_str = {'aeh', 'kon', 'sum'}; % Aeh Kon Sum tasks

% Find cue positions
perm_vector = cell(1, length(regr_cue_names)); % pre-allocate
for m_ind = 1:length(regr_cue_names)
    pos = length([perm_vector{:}]);
    perm_vector{m_ind} = [find(filter_condition(regr_cue_names{m_ind}, task_str{1})) ...
        find(filter_condition(regr_cue_names{m_ind}, task_str{2})) ...
        find(filter_condition(regr_cue_names{m_ind}, task_str{3}))] + pos;
end
perm_vector = [perm_vector{:}]; % convert to array

% Add model number to the regressor names
for l = 1:length(regr_cue_names)
    regr_cue_names{l} = strcat(regr_cue_names{l}, {['_' modelnames{l}]});
end

% Store all regressor names in the same array
regnames_all  = [regr_cue_names{:}];

% Permute array
regr_cue_names = regnames_all(perm_vector);

%% Check

assert(length(all_ts) == length(prefweight), ...
    'ts_data (all_ts) is meant to have one parameter set per 3dim "voxel" layer like prefweight, but all_ts = %i =/= prefweight = %i', ...
    length(all_ts), length(prefweight));

dispv(1, '    Generated data has %i voxel layers (possibly %i unique parameter sets)', ...
    length(prefweight), length(prefweight));

%% Create voxels

all_data = cell(1, length(all_ts)); % pre-allocate

for ts_ind = 1:length(all_ts)
    
    % Unpack current timeseries parameters
    ts_data = all_ts{ts_ind}; 
    
    % Add prespecified preferrred weights of trial and rest phases
    curr_vox_wmat = prefweight{ts_ind};
    
    % Add additional noise if specified
    curr_add_noise = 0;
    if exist('additionalnoise', 'var') && ~isempty(additionalnoise) && all(additionalnoise ~= 0)
        curr_add_noise = additionalnoise(ts_ind);
    end
    
    % Create voxels that respond only to a certain model
    voxel_data = struct();
    for run_ind = 1:length(ts_data) % runwise
        
        % Re-order runwise according to permutation vector
        ts_data(run_ind).X          = ts_data(run_ind).X(:, perm_vector);
        ts_data(run_ind).regr_names = regr_cue_names;
        
        % Initialise matrix in which we store the voxels
        voxel_data(run_ind).X          = []; 
        voxel_data(run_ind).voxelnames = [];
        voxel_data(run_ind).weights    = [];
        
        % Initialise row of current weights matrix)
        voxel_n = 0; row_cnt = 0;
        
        if cfg.saa.data_gen.partialresponsevoxels == 1 % NEW
            
            for m_ind = 1:length(models)
                
                modelname = sprintf('_%s', modelnames{m_ind});
                
                for c_ind = 1:length(task_str) % number of task types (classes)
                    
                    % Counters
                    row_cnt     = row_cnt + 1;
                    voxel_n     = voxel_n + 1;
                    
                    % Current variables
                    curr_cond   = task_str{c_ind};
                    curr_weight = curr_vox_wmat(row_cnt, c_ind);
                    type_filter = filter_condition(regr_cue_names, modelname, curr_cond);
                    
                    % Weigh and merge
                    voxels_w_merge = sum((ts_data(run_ind).X(:, type_filter) * curr_weight), 2);
                    
                    % Add noise if selected
                    curr_voxel = voxels_w_merge + curr_add_noise * randn(size(voxels_w_merge));
                    
                    % Add voxel to struct
                    voxel_data(run_ind).X(:, voxel_n)   = curr_voxel;
                    voxel_data(run_ind).voxelnames      = [voxel_data(run_ind).voxelnames {[curr_cond modelname]}];
                    voxel_data(run_ind).weights(row_cnt, c_ind) = curr_weight;
                    
                end % conditions (tasks)
            end % models
            
        else % every voxel responds to every task type (OLD VERSION)
            
            for m_ind = 1:length(models)
                
                modelname = sprintf('_%s', modelnames{m_ind});
                
                for cond_ind = 1:length(task_str) % number of tasks
                    
                    curr_cond   = task_str{cond_ind}; % Current task (aeh, kon, sum)
                    row_cnt  = row_cnt + 1;
                    voxel_n     = voxel_n + 1;
                    
                    aehw = curr_vox_wmat(row_cnt, 1);
                    konw = curr_vox_wmat(row_cnt, 2);
                    sumw = curr_vox_wmat(row_cnt, 3);
                    
                    filter_aeh = filter_condition(regr_cue_names, modelname, 'aeh');
                    filter_kon = filter_condition(regr_cue_names, modelname, 'kon');
                    filter_sum = filter_condition(regr_cue_names, modelname, 'sum');
                    
                    % Filter and weights
                    voxels_weighted_merged_aeh = sum((ts_data(run_ind).X(:, filter_aeh) * aehw), 2);
                    voxels_weighted_merged_kon = sum((ts_data(run_ind).X(:, filter_kon) * konw), 2);
                    voxels_weighted_merged_sum = sum((ts_data(run_ind).X(:, filter_sum) * sumw), 2);
                    
                    % Overall voxel signal (this copies the same timeseries
                    % for ever phase voxel regardless of class!!)
                    curr_voxel = voxels_weighted_merged_aeh + voxels_weighted_merged_kon + voxels_weighted_merged_sum;
                    
                    % Add noise if selected
                    curr_voxel = curr_voxel + curr_add_noise * randn(size(curr_voxel));
                    
                    % Add voxel to struct
                    voxel_data(run_ind).X(:, voxel_n)          = curr_voxel;
                    voxel_data(run_ind).voxelnames             = [voxel_data(run_ind).voxelnames {[curr_cond modelname]}];
                    voxel_data(run_ind).weights(row_cnt, :) = [aehw konw sumw];
                    
                end % conditions (tasks)
            end % models
        end % full or selective voxel response
    end % runs
    
    % Store current data
    all_data{ts_ind} = voxel_data;
    
end % layers
