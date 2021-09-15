% cfg_merged = afx_merge_complex_simple(all_similarity_cfg)
%
% Script to average corresponding betas in complex and simple conditions,
% e.g. AehSTrial10 and AehCTrial10.
%
% IN: all_similarity_cfg as a collection of cfgs (with results contained
%   inside) of the subject correlational analyses (see
%   run_timeresolved_analyses).
%
% OUT: 
%  'all_cond_cs' (135x135x54), storing averaged (via Fisher transformation)
%  correlation coefficients for all participants. Packed in cfg_merged.
% 
% Kai and Leonardo, 2019

function cfg_merged = afx_merge_complex_simple(all_similarity_cfg)  % all_conds_cs

subs_todo = all_similarity_cfg{1}.subs_todo;

for i = 1:length(all_similarity_cfg)
    subj = subs_todo(i);
    cond_merged_all{subj} = all_similarity_cfg{i}.saa.similarity_analysis_results.corr_results.cond_merged; %#ok<AGROW>
end

%% Set up variables
cfg1 = all_similarity_cfg{1}; % Take info from the first subject (it's the same for all participants)
con_unique_odd   = cfg1.saa.similarity_analysis_results.condition_unique.odd;
con_unique_even  = cfg1.saa.similarity_analysis_results.condition_unique.even;

%% Print some information in the command window

% label complex/simple pairs with same numbers
disp('Averaging complex and simple conditions')

% Init filter
unique_c_s_filter_odd = nan(size(con_unique_odd));

for u_ind = 1:length(con_unique_odd)
    curr_d = con_unique_odd{u_ind};
    c_s = curr_d(4); % rely on C/S is at position 4
    if strcmp(c_s, 'C')
        % find corresponding S condition
        s_string = [curr_d(1:3) 'S' curr_d(5:end)];
        s_ind = find(strcmp(con_unique_odd, s_string));
        if isempty(s_ind)
            error('Could not find simple condition %s, please check', s_string)
        elseif length(s_ind) > 2
            error('Found multiple instances for %s, please check', s_string)
        end
        % mark both with current index
        unique_c_s_filter_odd([u_ind s_ind]) = u_ind;
    elseif strcmp(c_s, 'S')
        % nothing to do
    else
        error('expecting C or S on 4th position, please check')
    end
end

% label complex/simple pairs with same numbers
unique_c_s_filter_even = nan(size(con_unique_even));

for u_ind = 1:length(con_unique_even)
    curr_d = con_unique_even{u_ind};
    c_s = curr_d(4); % rely on C/S is at position 4
    if strcmp(c_s, 'C')
        % find corresponding S condition
        s_string = [curr_d(1:3) 'S' curr_d(5:end)];
        s_ind = find(strcmp(con_unique_even, s_string));
        if isempty(s_ind)
            error('Could not find simple condition %s, please check', s_string)
        elseif length(s_ind) > 2
            error('Found multiple instances for %s, please check', s_string)
        end
        % mark both with current index
        unique_c_s_filter_even([u_ind s_ind]) = u_ind;
    elseif strcmp(c_s, 'S')
        % nothing to do
    else
        error('expecting C or S on 4th position, please check')
    end
end

%% Create data matrix with C/S conditions merged

fprintf('Calculating matrix with complex and simple conditions averaged within participants...'); layout_line_break(1);
all_cond_cs = cell(size(cond_merged_all));
for sb = subs_todo
    [cond_cs_z, all_cond_cs_labelnames_odd_orig, all_cond_cs_labelnames_even_orig ] = average_mat(atanh(cond_merged_all{sb}), ...
        unique_c_s_filter_odd, con_unique_odd, unique_c_s_filter_even, con_unique_even, @nanmean);
    all_cond_cs{sb} = tanh(cond_cs_z);
end
all_cond_cs = cat(3, all_cond_cs{:});

%% Get number of voxels

if strcmp(all_similarity_cfg{1}.data_type, 'simulated')
    saa_cfg = all_similarity_cfg{1}.saa.data_gen;
    if ~isfield(all_similarity_cfg{1}.saa.data_gen, 'total_n_voxels')
        for v_ind = saa_cfg.layer_number % 1:n
            try
                all_n_voxels(v_ind, 1) = length(saa_cfg.taskweights{v_ind}) * length(saa_cfg.phaseweights{v_ind}) + saa_cfg.n_mixedvoxels(v_ind); %#ok<AGROW>
            catch
                all_n_voxels(v_ind, 1) = length(saa_cfg.taskweights{1}) * length(saa_cfg.phaseweights{1}) + saa_cfg.n_mixedvoxels(1); %#ok<AGROW>
            end
        end
    else
        all_n_voxels = all_similarity_cfg{1}.saa.data_gen.total_n_voxels; % Cannot vary between subs
    end
    
    % Find bounds
    [vmin, vmax] = bounds(all_n_voxels);
    if vmin == vmax
        voxelrange = vmin;
    else
        voxelrange = [vmin vmax];
    end
    mean_n_voxels = mean(all_n_voxels);
    
else
    
    all_n_voxels = [];
    for i = 1:length(all_similarity_cfg)
        curr_sim_cfg = all_similarity_cfg{i};
        % Get number of voxels for current participant
        all_n_voxels(i) = curr_sim_cfg.saa.similarity_analysis_results.n_voxels; %#ok<AGROW>
    end
    
    % Get min and max number of voxels
    max_n_voxels    = max(all_n_voxels);
    min_n_voxels    = min(all_n_voxels);
    mean_n_voxels   = mean(all_n_voxels);
    voxelrange      = [min_n_voxels max_n_voxels];
end

%% Store result in cfg_step_2 struct

cfg_merged.general_info                              = cfg1;
cfg_merged.general_info.all_n_voxels                 = all_n_voxels;
cfg_merged.general_info.voxelrange                   = voxelrange;
cfg_merged.general_info.mean_n_voxels                = mean_n_voxels;
cfg_merged.group_data.merged_by_condition.dat        = cond_merged_all;
cfg_merged.group_data.merged_by_condition.descrip    = sprintf('Matrix for each participant containing the correlation coefficients\naveraged by experimental condition via Fisher transformation.\nNB: Complex and simple trials separate.');
cfg_merged.group_data.merged_by_condition_cs.dat     = all_cond_cs;
cfg_merged.group_data.merged_by_condition_cs.descrip = sprintf('Matrix for each participant containing the correlation coefficients\naveraged by experimental condition and complex/simple condition via Fisher transformation.');
cfg_merged.group_data.labelnames_cs_averaged.odd     = all_cond_cs_labelnames_odd_orig;
cfg_merged.group_data.labelnames_cs_averaged.even    = all_cond_cs_labelnames_even_orig;

%% Output end of script
% layout_line_break(1)
% layout_line('-', rep_symbol); layout_line_break(2);
disp('Simple and complex trial averaged')

end

