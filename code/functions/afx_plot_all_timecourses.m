%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot correlation coefficients for the odd runs against the even runs (time resolved)
% Averaged across participants and c/s merged
%
% Leonardo, 2019
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Disabled saving .fig for now, takes time and .fig are not really used,
% instead saves to svg using save_plot.m - Ingmar 23-07-21

%% Load workspace
function cfg_sim = afx_plot_all_timecourses(cfg_sim)

set(0, 'DefaultFigureVisible', 'off')
set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(0, 'DefaultAxesTickLabelInterpreter', 'none')

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
data_type                        = cfg_sim.general_info.data_type;
voxelrange                       = cfg_sim.general_info.voxelrange;
mean_n_voxels                    = cfg_sim.general_info.mean_n_voxels;
% %% Print some information in the command window
% indent                           = cfg.general_info.layout.header.indent;
% symbol                           = cfg.general_info.layout.header.symbol;
% rep_symbol                       = cfg.general_info.layout.header.rep_symbol;
% line_break                       = cfg.general_info.layout.header.line_break;

% layout_print_header('III. PLOT MEANS AND CONFIDENCE INTERVALS', symbol, rep_symbol, line_break)

if save_figures == 2
    fprintf('Figures will be saved in %s\n', targetdir);
else
    fprintf('Figures will not be saved\n');
end
layout_line_break(2);

%% Plot C/S merged block horizontal/vertical (i.e. the time course, odd vs even)

% Correlation to plot 2 options:
% 1) Single trial/resting state beta image in the middle of the block
% (plot_single_beta=1)
% 2) Averaged trials/resting state images (plot_single_beta=0)

% Select one trial beta and one resting beta in each block, odd runs

nameplot = 'time_resolved_correlation_all_time_series_complex_simple_averaged';
set_figure(cfg_sim.general_info, nameplot, [33   260   829   779]);
title(get(gcf, 'name'))

%% Source mat, source labels, target combination

% cnames    = x_labelnames;
% cnames_cs = x_labelnames_cs; % y and x names are the same

% filter_trial = filter_condition(cnames, 'Trial');
% filter_rest  = filter_condition(cnames, 'Rest');
% 
% filter_trial_cs = filter_condition(cnames_cs, 'Trial');
% filter_rest_cs  = filter_condition(cnames_cs, 'Rest');

% Source data
source_data    = all_means_cs;
source_labels  = x_labelnames_cs;

%% Select congruent conditions (Aeh-Aeh, Kon-Kon, Sum-Sum) for Trial and Resting state

% congruent_trial = select_congruent_conditions(source_data, source_labels, 'Trial');
% congruent_rest  = select_congruent_conditions(source_data, source_labels, 'Rest');
% 
% %% Select incongruent conditions (Aeh-Sum, Aeh-Kon, Kon-Sum, etc.) for Trial and Resting state
% 
% incongruent_trial = select_incongruent_conditions(source_data, source_labels, 'Trial');
% incongruent_rest  = select_incongruent_conditions(source_data, source_labels, 'Rest');

%% Adapt source data

source_data(source_data > 1-eps & source_data < 1+eps)=1-eps;
source_data(source_data > -1-eps & source_data < -1+eps)=-1+eps;

%% Plot curves

conds  = {'Aeh', 'Kon', 'Sum'};
conds2 = {'Trial', 'Rest'};
% PLOT SINGLE SURVES
subplot(2,2, [3:4])
colours = linspecer(18);
% colours = colours * 0.8; % make them a bit darker

hold all
curve_n = 0;
all_average_curves = struct();
for r = 1:length(conds)
    for c = 1:length(conds)
        for c2 = 1:length(conds2)
            curve_n = curve_n + 1;
            curr_cond_row  = conds{c};
            curr_state_row = conds2{c2};
            curr_cond_col  = conds{r};
            frow  =  filter_condition(source_labels, curr_cond_row, curr_state_row);
            fcol  =  filter_condition(source_labels, curr_cond_col);
            fmat = select_values_from_matrix(frow, fcol); % first filter: rows, second filter: columns
            current_filtered_mat = filter_matrix_with_matrix(source_data, fmat);
            
            % Average mat
            average_curve = tanh(nanmean(atanh(current_filtered_mat)));
            average_curve = average_curve(~isnan(average_curve));
            curve_name = [x_name '-' curr_cond_col '(Trial and Rest) X ' y_name '-' curr_cond_row '(' curr_state_row ')' ];
            
            % Calculate CIs
            mat = current_filtered_mat(:, ~all(isnan(current_filtered_mat)));
            mat = mat(any(~isnan(mat), 2), :);
            zmean = nanmean(atanh(mat));
            n = size(mat, 1);% get number of non-nan rows (n of curves to average)
            zCI95 = 1.96 * ( nanstd( atanh(mat), 0, 1) ./ sqrt(n));         
            lowerCI = tanh(zmean-zCI95);
            upperCI = tanh(zmean+zCI95);
             
            % Store results
            all_average_curves(curve_n).curve               = average_curve;
            all_average_curves(curve_n).lowerCI             = lowerCI;
            all_average_curves(curve_n).upperCI             = upperCI;
            all_average_curves(curve_n).curve_name          = curve_name;
            all_average_curves(curve_n).curve_number        = curve_n;
            all_average_curves(curve_n).curve_name_for_plot = [num2str(curve_n) '-' curve_name];
            all_average_curves(curve_n).mats_trace          = current_filtered_mat;
            all_average_curves(curve_n).mats_to_average     = mat;
            
            % Plot curve
            if sum(~isnan(average_curve))
                hold on
                h(curve_n) = plot(average_curve,'LineWidth', 2, 'Color', colours(curve_n, :));
                plot(lowerCI, '-.','LineWidth', 0.2, 'Color', colours(curve_n, :))
                plot(upperCI, '-.','LineWidth', 0.2, 'Color', colours(curve_n, :))
                hold off
            end
        end
    end
end


%% Set up plot 

% Time points 
xline(1, ':', 'cue', 'LabelVerticalAlignment', 'bottom');
xline(2, ':', 'trial onset', 'LabelVerticalAlignment', 'bottom');
xline(27, ':', 'trial end', 'LabelVerticalAlignment', 'bottom');
xline(28, ':', 'resting state', 'LabelVerticalAlignment', 'bottom');
xline(45, ':', 'end block', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
% Ref line
line(xlim(), [0,0], 'LineWidth', 0.2, 'Color', 'k');
xlim([1 45])

% Labels
gap = 5; % 1: no gap
set(gca, 'Xtick', 1:gap:length(source_labels), ...
    'XTickLabel', source_labels(1:gap:length(source_labels)), 'XTickLabelRotation', 90);
% Legend

if exist('h', 'var')
    legend(h, all_average_curves.curve_name_for_plot); % all_average_curves.curve_name;
end

set(legend, 'Location', 'eastoutside'); % position legend outside of the graph
set(legend, 'Interpreter', 'none'); % set interpreter to 'none' to avoid weird random characters being added
ylabel('mean(RHOs)');
xlabel(sprintf('volumes (%s runs, c/s averaged)', x_name));
% Add figure name (depending on line style too)
figname = fullfile(targetdir, sprintf('fig-%s', get(gcf, 'name')));

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
for cb = 1:length(all_average_curves)
    curr_mat = all_average_curves(cb).mats_trace;
    selected_blocks(~isnan(curr_mat)) = cb;
    % Get coordinates of the centre of the block (to write the curve number
    [row_idx, column_idx] = find(selected_blocks == cb);
    xtext = round(mean(row_idx));
    ytext = round(mean(column_idx));
    coords(cb, :) = [xtext ytext];
end
selected_blocks(contains(source_labels, 'Cue'), :) = nan; % Set cues on the rows to 0
selected_blocks(:, contains(source_labels, 'Cue')) = nan; % Set column on the rows to 0
imagesc(selected_blocks, 'AlphaData', ~isnan(selected_blocks))
% Plot curve number on matrix
for cb = 1:length(all_average_curves)
    text(coords(cb, 2), coords(cb, 1), num2str(cb)) % NB: text position parameters x and y correspond to columns and rows, respectively (i.e. swap them)
end

title(gca, 'Selected data')
colormap(gca, colours)
colorbar
axis equal
axis tight

% Labels
xlabel(sprintf('%s runs', x_name));
ylabel(sprintf('%s runs', y_name));

%% Add suptitle
if strcmp(data_type, 'wholebrain')|| strcmp(data_type, 'empirical')
    voxelinfo = sprintf('whole brain - n_voxels[min=%i,max=%i,mean=%g]', voxelrange(1), voxelrange(2), round(mean_n_voxels));
elseif strcmp(data_type, 'ROI')
    voxelinfo = sprintf('ROI - n_voxels[min=%i,max=%i,mean=%g]', voxelrange(1), voxelrange(2), round(mean_n_voxels));
elseif strcmp(data_type, 'simulated')
    voxelinfo = sprintf('simulated - n_voxels[%i]', voxelrange(1));
end

suptitle(voxelinfo);

%% Save figure if selected
disp(['Storing: ' figname]);
if save_figures == 2
    save_plot(figname, [], 0, 0, 1)
    %save_fig(figname);
end

%% Store useful info in the cfg structure
cfg_sim.timecourse.all_average_curves = all_average_curves;

end




