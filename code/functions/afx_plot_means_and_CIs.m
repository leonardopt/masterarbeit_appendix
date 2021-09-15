%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data Analysis and Visualisation: calculate and plot odd means and confidence intervals (time-resolved respectively to even)
%
% First Part (keeps simple and complex separate)
% 1. Calculate average correlation matrices across participants via Fisher transformation
%    - output variable --> 'all_means'
%    - plot 1, display means using image with scaled colours
% 2. Calculate sandard error of the mean, standard deviation, and 95% confidence intervals
%    - plot 2, display CIs using image with scaled colours
%    - plot 3, scatterplot of CI95s against means (mean+/- CI and mean)

% Second Part: do exactly the same, but this time using the C/S averaged data
% 	- NB: needs to load variable 'all_cond_cs' from 'all_cond_cs.mat' to work


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cfg_sim = afx_plot_means_and_CIs(cfg_sim)

set(0, 'DefaultFigureVisible', 'off')

save_figures                     = cfg_sim.general_info.save_figures;
targetdir                        = cfg_sim.general_info.target_directory;
cond_merged_all                  = cfg_sim.group_data.merged_by_condition.dat;
all_cond_cs                      = cfg_sim.group_data.merged_by_condition_cs.dat;
odd_labelnames                   = cfg_sim.general_info.saa.similarity_analysis_results.labelnames.odd;
even_labelnames                  = cfg_sim.general_info.saa.similarity_analysis_results.labelnames.even;
odd_labelnames_cs                = cfg_sim.group_data.labelnames_cs_averaged.odd ;
even_labelnames_cs               = cfg_sim.group_data.labelnames_cs_averaged.even;
voxelrange                       = cfg_sim.general_info.voxelrange;
mean_n_voxels                    = cfg_sim.general_info.mean_n_voxels;
data_type                        = cfg_sim.general_info.data_type;


%% Print some information in the command window

if save_figures == 2
    fprintf('Figures will be saved in %s\n', targetdir);
else
    fprintf('Figures will not be saved\n');
end
layout_line_break(2);



%% %%%%%%%%%%%%%%% COMPLEX AND SIMPLE SEPARATE %%%%%%%%%%%%%%% %%

fprintf('1. Analysing correlations - complex and simple conditions kept separate');layout_line_break(2);
layout_indent(10); disp('Calculating means and CI95s...'); layout_line_break(1);

%% Average correlation matrices across participants

mean_mat_name ='similarity_matrix_mean_rho_complex_simple_separate';
CI_mat_name   =  'similarity_matrix_CI95_complex_simple_separate';
[all_means, all_stds, all_sems, all_CI95, odd_labelnames_polished, even_labelnames_polished] = afx_tool_plot_similarity_matrices(cfg_sim, mean_mat_name, CI_mat_name, cond_merged_all, odd_labelnames, even_labelnames);


%% Plot mean vs CI95 and vs mean+CI95

layout_indent(10); disp('Plotting mean vs CI95...'); layout_line_break(1);

% Create mask for colors
odd_labelnames = odd_labelnames' ;
even_labelnames = even_labelnames';

mean_vs_CI_name ='mean_vs_CI95_complex_simple_separate';
afx_tool_plot_mean_vs_CI95(cfg_sim, mean_vs_CI_name, all_means, all_CI95, odd_labelnames, even_labelnames)


%% %%%%%%%%%%%%%%% COMPLEX AND SIMPLE AVERAGED %%%%%%%%%%%%%%% %%

% C/S merged ------ Average correlation matrices across participants
fprintf('2. Analysing correlations - complex and simple conditions averaged'); layout_line_break(2);
layout_indent(10); disp('Calculating means and CI95s...'); layout_line_break(1);

%% Average correlation matrices across participants

mean_mat_name_cs ='similarity_matrix_mean_rho_complex_simple_averaged';
CI_mat_name_cs   = 'similarity_matrix_CI95_complex_simple_averaged';
[all_means_cs, all_stds_cs, all_sems_cs, all_CI95_cs,odd_labelnames_cs_polished, even_labelnames_cs_polished] = afx_tool_plot_similarity_matrices(cfg_sim, mean_mat_name_cs, CI_mat_name_cs, all_cond_cs, odd_labelnames_cs, even_labelnames_cs);



%% C/S merged ------ Plot mean vs CI95 and vs mean+CI95

layout_indent(10); disp('Plotting mean vs CI95...'); layout_line_break(1);

% Create mask for colors
odd_labelnames_cs = odd_labelnames_cs' ;
even_labelnames_cs = even_labelnames_cs';


mean_vs_CI_name_cs ='mean_vs_CI95_complex_simple_averaged';
afx_tool_plot_mean_vs_CI95(cfg_sim, mean_vs_CI_name_cs, all_means_cs, all_CI95_cs, odd_labelnames_cs, even_labelnames_cs);




%% Update cfg structa
cfg_sim.general_info.settings.plotting.x_name  = 'even';
cfg_sim.general_info.settings.plotting.y_name  = 'odd';
cfg_sim.general_info.settings.plotting.x_labelnames_polished    = even_labelnames_polished;
cfg_sim.general_info.settings.plotting.y_labelnames_polished    = odd_labelnames_polished;
cfg_sim.general_info.settings.plotting.x_labelnames_cs_polished = even_labelnames_cs_polished;
cfg_sim.general_info.settings.plotting.y_labelnames_cs_polished = odd_labelnames_cs_polished;

cfg_sim.second_level_statistics.merged_by_condition.all_means    = all_means;
cfg_sim.second_level_statistics.merged_by_condition.all_stds     = all_stds;
cfg_sim.second_level_statistics.merged_by_condition.all_sems     = all_sems;
cfg_sim.second_level_statistics.merged_by_condition_cs.all_means = all_means_cs;
cfg_sim.second_level_statistics.merged_by_condition_cs.all_stds  = all_stds_cs;
cfg_sim.second_level_statistics.merged_by_condition_cs.all_sems  = all_sems_cs;
cfg_sim.time_log.step_3                                          = datestr(now);


%% Output end of script
layout_line_break(1)
layout_line('=', 80); layout_line_break(2)

end

