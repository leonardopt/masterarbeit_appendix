% [all_means, all_stds, all_sems, all_CI95, odd_labelnames_polished, ...
%     even_labelnames_polished] = afx_tool_plot_similarity_matrices(cfg_sim, ...
%     mean_mat_name, CI_mat_name, data, odd_labelnames, even_labelnames)
% 
% Plotting function by Leonardo, 2019
% 
% Disabled saving .fig for now, takes time and .fig are not really used,
% instead saves to svg using save_plot.m - Ingmar 23-07-21

function [all_means, all_stds, all_sems, all_CI95, odd_labelnames_polished, ...
    even_labelnames_polished] = afx_tool_plot_similarity_matrices(cfg_sim, ...
    mean_mat_name, CI_mat_name, data, odd_labelnames, even_labelnames)

set(0, 'DefaultFigureVisible', 'off')

save_figures                     = cfg_sim.general_info.save_figures;
targetdir                        = cfg_sim.general_info.target_directory;
voxelrange                       = cfg_sim.general_info.voxelrange;
mean_n_voxels                    = cfg_sim.general_info.mean_n_voxels;
data_type                        = cfg_sim.general_info.data_type;

% Create string with voxel number info for the titles
if ~strcmp(data_type, 'simulated')
    voxelinfo = sprintf('n_voxels[min=%i,max=%i,mean=%g]', voxelrange(1), voxelrange(2), round(mean_n_voxels));
end

%% Average correlation matrices across participants
% Get mean according to type of input data
% 1) Perform Fisher z transformation: corr coefficients r --> z
% 2) Calculate mean in terms of z
% 3) Transform back z --> r

if iscell(data)
    
    % First eliminiate empty cells from cond_merged_all, if there are any
    data = data(~cellfun('isempty', data));
    all_data = cell2mat(reshape(data, [1,1, size(data, 2)]));
    all_means = tanh((nanmean(atanh(all_data), 3)));
    
elseif isa(data, 'double')
    
    all_means = tanh((nanmean(atanh(data), 3)));
    all_data = data;
end

%% Show result

set_figure(cfg_sim.general_info, mean_mat_name, [103   186   972   853]);

imagesc(all_means)
mark_nans(all_means)


if strcmp(data_type, 'simulated')
    name_fig = get(gcf, 'name');
    title( sprintf('%s\n[%s] voxels', name_fig, num2str(voxelrange(1))));
else
    name_fig = get(gcf, 'name');  
    title(sprintf('%s\n%s', name_fig, voxelinfo));
end

colorbar
% Set labels
[odd_labelnames_polished, even_labelnames_polished] = set_labels_for_matrices(odd_labelnames, even_labelnames);

% Print in command window
layout_indent(10); fprintf('Mean via Fisher transformation -- done'); layout_line_break(1);

% Save figures, if selected by user
figname = fullfile(targetdir, sprintf('fig-%s', get(gcf, 'name')));
if save_figures == 2
    disp(['Storing: ' figname]);
    save_plot(figname, [], 0, 0, 1)
    % save_fig(figname); % without TDT: save_fig, but then only .fig (check help for .eps and .png)
end
layout_line_break(1);

%% Calculate Standard Error of the Mean (SEM), standard deviation (STD) and 95% Confidence Intervals (CI95)

all_counts = size(all_data, 3);
all_stds = tanh(nanstd(atanh(all_data), 0, 3)); % uses N-1 for std
all_sems =  tanh( nanstd( atanh(all_data), 0, 3) ./ sqrt(all_counts) ); % uses N-1 for std
all_CI95 =  tanh( nanstd( 1.96 * atanh(all_data), 0, 3) ./ sqrt(all_counts) ); % uses N-1 for std

set_figure(cfg_sim.general_info, CI_mat_name,[103   186   972   853]);

imagesc(all_CI95);
mark_nans(all_CI95)

if strcmp(data_type, 'simulated')
    name_fig = get(gcf, 'name'); 
    title( sprintf('%s\n[%s] voxels', name_fig, num2str(voxelrange(1))));
    
elseif strcmp(data_type, 'empirical') || strcmp(data_type, 'wholebrain')
    name_fig = get(gcf, 'name');  title(sprintf('%s\n%s', name_fig, voxelinfo));
    title( sprintf('%s\nwhole brain - [%s] voxels', name_fig, num2str(voxelrange(1))));

elseif strcmp(data_type, 'ROI')
      name_fig = get(gcf, 'name');  title(sprintf('%s\n%s', name_fig, voxelinfo));
    title( sprintf('%s\nROI - [%s] voxels', name_fig, num2str(voxelrange(1))));

end

colorbar;
set_labels_for_matrices(odd_labelnames, even_labelnames);

% Print in command window
layout_indent(10); fprintf('CI95 via Fisher transformation -- done'); layout_line_break(2);

% Save figures, if selected by user
figname = fullfile(targetdir, sprintf('fig-%s', get(gcf, 'name')));
if save_figures == 2
    disp(['Storing: ' figname]);
    save_plot(figname, [], 0, 0, 1)
    % save_fig(figname); % without TDT: save_fig, but then only .fig (check help for .eps and .png)
end
layout_line_break(1);

end 