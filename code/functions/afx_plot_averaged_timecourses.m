%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot correlation coefficients for the odd runs against the even runs (time resolved)
% Averaged across participants and c/s merged
% Average trial/resting state data from the same odd/even condition and from different odd/even conditions
%
% 4 combinations:
% - averaged trials same conditions odd and even (3 in total)
%       - e.g. mean(AehTrials odd x AehTrials even, KonTrials odd x KonTrials even, SumTrials odd x SumTrials even)
% - averaged trials different conditions (odd) vs even (6 in total)
%       - e.g. mean(AehTrials odd x KonTrials even, KonTrials odd x SumTrials even, SumTrials odd x AehTrials even)
% - averaged resting states same conditions (odd) vs even (3 in total)
% - averaged resting states different conditions (odd) vs even (6 in total)
%
% Two options:
% 1) Trials to average are single trials taken from the blocks
%   - e.g. mean(AehTrial26 odd x AehTrials even, KonTrial26 odd x KonTrials even, SumTrial26 odd x SumTrials even)
%   - plot source data, traces ('horizontal'), graph with 4 curves (one per combination)
%   WARNING: check if it works - Leonardo 21/08/2019
%
%
% 2) Trials to average are all trials taken from the blocks (see above)
%   - plot graph with 4 curves (one per combination)
%
%  TO ADD: graph with source data and traces in cond 2 - Leonardo 21/08/2019
% use filter computed in 3.1.?
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Disabled saving .fig for now, takes time and .fig are not really used,
% instead saves to svg using save_plot.m - Ingmar 23-07-21

function cfg_sim = afx_plot_averaged_timecourses(cfg_sim)

set(0, 'DefaultLegendInterpreter', 'none')
set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultFigureVisible', 'off')

%% Set variables

targetdir                        = cfg_sim.general_info.target_directory;
save_figures                     = cfg_sim.general_info.save_figures;
x_name                           = cfg_sim.general_info.settings.plotting.x_name;
y_name                           = cfg_sim.general_info.settings.plotting.y_name;
% x_labelnames                     = cfg_sim.general_info.settings.plotting.x_labelnames_polished;
% y_labelnames                     = cfg_sim.general_info.settings.plotting.y_labelnames_polished;
x_labelnames_cs                  = cfg_sim.general_info.settings.plotting.x_labelnames_cs_polished;
% y_labelnames_cs                  = cfg_sim.general_info.settings.plotting.y_labelnames_cs_polished;
% all_means                        = cfg_sim.second_level_statistics.merged_by_condition.all_means;
all_means_cs                     = cfg_sim.second_level_statistics.merged_by_condition_cs.all_means;
all_means_cs_allsubs             = cfg_sim.group_data.merged_by_condition_cs.dat;
data_type                        = cfg_sim.general_info.data_type;
voxelrange                       = cfg_sim.general_info.voxelrange;
mean_n_voxels                    = cfg_sim.general_info.mean_n_voxels;

%% Filter conditions 

% cnames    = x_labelnames;
% cnames_cs = x_labelnames_cs; % y and x names are the same
% 
% filter_trial = filter_condition(cnames, 'Trial');
% filter_rest  = filter_condition(cnames, 'Rest');
% 
% filter_trial_cs = filter_condition(cnames_cs, 'Trial');
% filter_rest_cs  = filter_condition(cnames_cs, 'Rest');

% Source data
source_data    = all_means_cs;
source_labels  = x_labelnames_cs;

%% Select congruent conditions (Aeh-Aeh, Kon-Kon, Sum-Sum) for Trial and Resting state

for i  = 1:size(all_means_cs_allsubs, 3)
    congruent_trial{i} = select_congruent_conditions(all_means_cs_allsubs(:,:,i), source_labels, 'Trial');
    congruent_rest{i}  = select_congruent_conditions(all_means_cs_allsubs(:,:,i), source_labels, 'Rest');
end

% Select congruent for diagonal plotting
% from all participants
congruent_trial_and_rest = {};
for i  = 1:size(all_means_cs_allsubs, 3)
    congruent_trial_and_rest{i} = select_congruent_conditions(all_means_cs_allsubs(:,:,i), source_labels, 'Trial', 'Rest');
end

%% Select incongruent conditions (Aeh-Sum, Aeh-Kon, Kon-Sum, etc.) for Trial and Resting state

for i  = 1:size(all_means_cs_allsubs, 3)
    incongruent_trial{i} = select_incongruent_conditions(all_means_cs_allsubs(:,:,i), source_labels, 'Trial');
    incongruent_rest{i}  = select_incongruent_conditions(all_means_cs_allsubs(:,:,i), source_labels, 'Rest');
end

% Select incongruent for diagonal plotting from all participants 
incongruent_trial_and_rest = {};
for i  = 1:size(all_means_cs_allsubs, 3)
    incongruent_trial_and_rest{i} = select_incongruent_conditions(all_means_cs_allsubs(:,:,i), source_labels, 'Trial', 'Rest');
end


%% Plot C/S merged block horizontal/vertical (i.e. the time course, odd vs even)

% Set figure
nameplot_h = 'time_resolved_correlation_horizontal_averaged_time_series_complex_simple_averaged';
set_figure(cfg_sim.general_info, nameplot_h, [33   260   829   779]);
 
%% Average and plot
cond_n = 3; % number of conditions (6 when complex and simple are not averaged)
colors = [  0.9296    0.3098    0.3098
            0.9949    0.5666    0.5666
            0.3059    0.3059    0.9105
            0.5140    0.5140    0.9916];

symbols = {'-', '-', '-','-'};

mat_to_average = {congruent_trial, incongruent_trial, congruent_rest, incongruent_rest};
curv_names     = {'congruent_trial', 'incongruent_trial', 'congruent_rest', 'incongruent_rest'};

subplot(2, 2, [3:4])
for cv = 1:length(mat_to_average)
  
   [average_curves{cv}, lowerCI{cv}, upperCI{cv}] = average_selected_blocks(mat_to_average{cv}, cond_n);
    xaxis_values = 1:size(average_curves{cv}, 2);
    %     y1 = lowerCI{cv};
    %     y2 = upperCI{cv};
    %
    %     x_plot = [x, fliplr(x)];
    %     y_plo§t = [y1, flipud(y2)];
    %     fill(x_plot, y_plot, 1,'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.4);
    hold on
    %     plot(average_curves{cv}, symbols{cv}, 'color', colors(cv, :))
    %     plot(lowerCI{cv}, ':', 'color', colors(cv, :))
    %     plot(upperCI{cv}, ':', 'color', colors(cv, :))
    
    if sum(~isnan(average_curves{cv}))  
        h(cv) = plot(average_curves{cv}, symbols{cv}, 'LineWidth', 2, 'color', colors(cv, :));
        plot(lowerCI{cv}, '-.','LineWidth', 0.2, 'color', colors(cv, :))
        plot(upperCI{cv}, '-.','LineWidth', 0.2, 'color', colors(cv, :))
    end 
    %
    %     xaxis_fill = [xaxis_values, fliplr(xaxis_values)];
    %     yaxis_fill = [ upperCI{cv}, flipud(lowerCI{cv})];
    %     fill(xaxis_fill, yaxis_fill, [1 0.84 0.21], 'facecolor',  ... % lower sat yellow/orange
    %             [1 0.84 0.21], 'edgecolor', [1 0.84 0.21], 'facealpha', 0.15);
end

%% Plot settings

% Time points 
xline(1, ':', 'cue', 'LabelVerticalAlignment', 'bottom');
xline(2, ':', 'trial onset', 'LabelVerticalAlignment', 'bottom');
xline(27, ':', 'trial end', 'LabelVerticalAlignment', 'bottom');
xline(28, ':', 'resting state', 'LabelVerticalAlignment', 'bottom');
xline(45, ':', 'end block', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
% Ref line
line(xlim(), [0,0], 'LineWidth', 0.2, 'Color', 'k');
xlim([1 45])

% Legend
set(legend, 'Interpreter', 'none')
if exist('h', 'var')
    legend(h, curv_names)
end 

% Labels
xlabel(sprintf('%s runs timeseries', x_name));
ylabel(sprintf('%s runs - Mean(RHOs)', y_name));

t_start_even = find(~cellfun(@isempty, strfind(x_labelnames_cs, '00')));
gap = 10;
ticks_odd_av = unique([1:gap:length(x_labelnames_cs), t_start_even]);

set(gca, 'Xtick', ticks_odd_av, ...
    'XTickLabel', x_labelnames_cs(ticks_odd_av), 'XTickLabelRotation', 90,'FontSize',7);
% title({'Correlation cofficient averaged', 'by task congruency (odd runs) over time '});

%% Plot source data and selected blocks

% Source data
subplot(2,2,1)
imagesc(source_data)
title(gca, 'Source data')
colorbar
axis equal
axis tight

% Labels
xlabel(sprintf('%s runs', x_name));
ylabel(sprintf('%s runs', y_name));

% Selected blocks
subplot(2,2,2)
selected_blocks = nan(size(source_data));
for cb = 1:length(mat_to_average)
    curr_mat = mat_to_average{cb};curr_mat = curr_mat{1}; % take first (all the same)
    selected_blocks(~isnan(curr_mat)) = cb;
end
selected_blocks(contains(source_labels, 'Cue'), :) = nan; % Set cues on the rows to 0
selected_blocks(:, contains(source_labels, 'Cue')) = nan; % Set column on the rows to 0
imagesc(selected_blocks, 'AlphaData', ~isnan(selected_blocks))
title(gca, 'Selected data')
colormap(gca, colors)
colorbar
axis equal
axis tight

% Labels
xlabel(sprintf('%s runs', x_name));
ylabel(sprintf('%s runs', y_name));

% Colorbar
cb = colorbar('YTick',[0 0.5], 'YTicklabel',[], 'FontSize', 17, 'FontName', 'Calibri');
text(152, 70, sprintf('incongruent_rest\n\n\n\ncongruent_rest\n\n\n\nincongruent_trial\n\n\n\ncongruent_trial'), 'FontSize', 10, 'Rotation', 0);

figname = fullfile(targetdir, sprintf('fig-%s', get(gcf, 'name')));

%% Add suptitle

if strcmp(data_type, 'wholebrain') || strcmp(data_type, 'empirical')
    voxelinfo = sprintf('whole brain - n_voxels[min=%i,max=%i,mean=%g]', voxelrange(1), voxelrange(2), round(mean_n_voxels));
elseif strcmp(data_type, 'ROI')
    voxelinfo = sprintf('ROI - n_voxels[min=%i,max=%i,mean=%g]', voxelrange(1), voxelrange(2), round(mean_n_voxels));
elseif strcmp(data_type, 'simulated')
    voxelinfo = sprintf('simulated - n_voxels[%i]', voxelrange(1));
end

suptitle(voxelinfo);

%% Save figure if selected

if save_figures == 2
    save_plot(figname, [], 0, 0, 1)
    % save_fig(figname);
end

%% PLOT DIAGONAL

% Set figure
nameplot_d = 'time_resolved_correlation_diagonal_averaged_time_series_complex_simple_averaged';
set_figure(cfg_sim.general_info, nameplot_d, [33   260   829   779]);
subplot(2, 2, [3:4])

%% Average and plot diagonal

colors2 = [  0.9296         0.3098     0.3098
    0.3059    0.3059    0.9105];

symbols2 = {'-', '-'};

mat_to_average_diagonal = {congruent_trial_and_rest, incongruent_trial_and_rest};
curv_names     = {'congruent_diagonal', 'incongruent_diagonal'};

subplot(2, 2, [3:4])
for cv = 1:length(mat_to_average_diagonal)
    
    [average_curves2{cv}, lowerCI2{cv}, upperCI2{cv}] = average_selected_blocks_diagonal(mat_to_average_diagonal{cv}, cond_n);
    
    hold on
    if sum(~isnan(average_curves2{cv}))
        h(cv) = plot(average_curves2{cv}, symbols2{cv}, 'LineWidth', 2, 'color', colors2(cv, :));
        plot(lowerCI2{cv}, '-.','LineWidth', 0.2, 'color', colors2(cv, :))
        plot(upperCI2{cv}, '-.','LineWidth', 0.2, 'color', colors2(cv, :))
    end    
end

%% Plot settings 

% Time points 
xline(1, ':', 'cue', 'LabelVerticalAlignment', 'bottom');
xline(2, ':', 'trial onset', 'LabelVerticalAlignment', 'bottom');
xline(27, ':', 'trial end', 'LabelVerticalAlignment', 'bottom');
xline(28, ':', 'resting state', 'LabelVerticalAlignment', 'bottom');
xline(45, ':', 'end block', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
% Ref line
line(xlim(), [0,0], 'LineWidth', 0.2, 'Color', 'k');
xlim([1 45])

% Legend
set(legend, 'Interpreter', 'none')

if exist('h', 'var')
    legend(h, curv_names)
end

% Labels
xlabel(sprintf('%s runs timeseries', x_name));
ylabel(sprintf('%s runs - Mean(RHOs)', y_name));

t_start_even = find(~cellfun(@isempty, strfind(x_labelnames_cs, '00')));
gap = 10;
ticks_odd_av = unique([1:gap:length(x_labelnames_cs), t_start_even]);

set(gca, 'Xtick', ticks_odd_av, ...
    'XTickLabel', x_labelnames_cs(ticks_odd_av), 'XTickLabelRotation', 90,'FontSize',7);
% title({'Correlation cofficient averaged', 'by task congruency (odd runs) over time '});

%% Plot source data and selected blocks
% Source data
subplot(2,2,1)
imagesc(source_data)
title(gca, 'Source data')
colorbar
axis equal
axis tight
% Labels
xlabel(sprintf('%s runs', x_name));
ylabel(sprintf('%s runs', y_name));
% Selected blocks

subplot(2,2,2)

selected_blocks2 = reshape(1:numel(source_data), size(source_data));

% Index diagonal2
diag1 = diag(selected_blocks2);
diag2 = diag(selected_blocks2, 45);
diag3 = diag(selected_blocks2, -45);
diag4 = diag(selected_blocks2, 90)';
diag5 = diag(selected_blocks2, -90)';

% Plot traces on matrix 
selected_blocks2(diag1) = 0;
selected_blocks2(diag2) = 0;
selected_blocks2(diag3) = 0;
selected_blocks2(diag4) = 0;
selected_blocks2(diag5) = 0;

selected_blocks2(selected_blocks2~=0) = nan;

% Plot traces in different colours depending on whether they are congruent
% or incongruent conditions 
selected_blocks2(diag1) = 1;
selected_blocks2(diag2) = 2;
selected_blocks2(diag3) = 2;
selected_blocks2(diag4) = 2;
selected_blocks2(diag5) = 2;

imagesc(selected_blocks2, 'AlphaData', ~isnan(selected_blocks2))
title(gca, 'Selected data')
colormap(gca, colors2)
colorbar
axis equal
axis tight

% Labels
xlabel(sprintf('%s runs', x_name));
ylabel(sprintf('%s runs', y_name));

% Lines
l1 = yline(45, ':'); l1.Color=[0,0,0,0.2];
l2 = yline(90, ':'); l2.Color=[0,0,0,0.2];
l3 = xline(45, ':'); l3.Color=[0,0,0,0.2];
l4 = xline(90, ':'); l4.Color=[0,0,0,0.2];

% Colourbar
cb = colorbar('YTick',[0 0.5], 'YTicklabel',[], 'FontSize', 17, 'FontName', 'Calibri');
text(160, 130, 'congruent                     incongruent', 'FontSize', 10, 'Rotation', 90);

figname2 = fullfile(targetdir, sprintf('fig-%s', get(gcf, 'name')));

%% Add suptitle

if strcmp(data_type, 'wholebrain') || strcmp(data_type, 'empirical')
    voxelinfo = sprintf('whole brain - n_voxels[min=%i,max=%i,mean=%g]', voxelrange(1), voxelrange(2), round(mean_n_voxels));
elseif strcmp(data_type, 'ROI')
    voxelinfo = sprintf('ROI - n_voxels[min=%i,max=%i,mean=%g]', voxelrange(1), voxelrange(2), round(mean_n_voxels));
elseif strcmp(data_type, 'simulated')
    voxelinfo = sprintf('simulated - n_voxels[%i]', voxelrange(1));
end

suptitle(voxelinfo);

%% Save figure if selected
if save_figures == 2
    save_plot(figname2, [], 0, 0, 1)
    %save_fig(figname2);
end

%% Store useful info in the cfg_sim structure
cfg_sim.timecourse.selected_blocks = selected_blocks;
cfg_sim.timecourse.source_data     = source_data;
cfg_sim.timecourse.source_labels   = source_labels;

end