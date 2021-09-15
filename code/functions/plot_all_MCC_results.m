% function [resfig] = plot_all_MCC_results(result, targetfile, datatype)
%
% Constructs subplots that show hypothesis testing, p-values, mean group
% accuracies and confidence intervals for all training x testing
% timepoints for group condition 1 vs group condition 2 t-test results.
%
% This function only plots. The calculation is done in
% afx_grouplevel_decoding_stats. This figure shows all relevant stats in 1
% single file for easy comparison.
%
% Shows all multiple multiple comparison correction options alongside each
% other for a given brain space (e.g. ROI or whole-brain):
%       0: uncorr
%       1: bonf.
%       2: FDR 0.001
%       3: signperm_test clustersize reconstructed from individual quandrant matrices
%       4: signperm_test clustermass reconstructed from individual quandrant matrices
%       5: TFCE
%       6: TFCE reconstructed from individual quandrant matrices
%
% OUT
%   resfig{do_MCC_pos}: handle to figure for test at position do_MCC (not for the do_MCC value)
%   targetfile_fig{do_MCC_pos}: file for resfig{do_MCC_pos}
%
% Ingmar, 08-09-21

function [resfig] = plot_all_MCC_results(result, targetfile, datatype)
%% Settings

% Empirical or simulated?
if ~exist('datatype', 'var') || isempty(datatype)
    datatype = 'empirical'; disp('Assuming real data');
end

% No group name
if ~isfield(result, 'group_resname')
    result.group_resname = 'Group-level decoding results';
end

% Timing info (ticks, labels)
[customticks, customlabels] = afx_set_timings(0);

% Colors
teal = [43 147 166]; teal = teal/255;% lightblue = [53/255 119/255 233/255];
orange = [255 156 62]; orange = orange/255;

task_rest_tp = [14 14]; % <---- SET VALUES FOR TRAIN TP SLICES

%% Plot multiple comparison corrections

% Figure size
n_rows = 7; n_cols = 4; n_subplots = n_rows * n_cols;
fig_sz = [265 -2519 1716 3367] *1.5; % [1000      -2500       1788        3508]; % A4 format
font_sz = 8;

% Create full figure
resfig = figure('name', result.group_resname, 'Position', fig_sz);

%% Make subplots 1

do_MCC = 1; % bonf
[subplot_stat, subplot_titles, ~, ~] = select_MCC_type(do_MCC, result, targetfile); % Result/MCC type

% Mean
sp_ind = 1;
curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
imagesc(result.mean(1:43, 1:43));
decoding_plot_settings(curr_sp_hdl, datatype, font_sz,  customticks, customlabels); % Standard settings
title(subplot_titles{1}, 'Interpreter', 'none');

% Alter color limits
if isfield(result, 'extend_caxis') && ~isempty(result.extend_caxis)
    if contains(subplot_titles(sp_ind), 'width')
        caxis(result.extend_caxis(2, :))
    else
        caxis(result.extend_caxis(1, :))
    end
end

% CI95 1
sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
imagesc(subplot_stat{2}(1:43, 1:43));
bar_h = decoding_plot_settings(curr_sp_hdl, datatype, font_sz,  customticks, customlabels); % Standard settings
title(subplot_titles{2}, 'Interpreter', 'none');
xlabel(bar_h, 'Accuracy CI95 (%)', 'FontSize', font_sz); % overwrite label

% H 1
sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
imagesc(subplot_stat{3}(1:43, 1:43));
bar_h = decoding_plot_settings(curr_sp_hdl, datatype, font_sz,  customticks, customlabels); % Standard settings
colormap(curr_sp_hdl, [teal; 1 1 1; orange]);
xlabel(bar_h, 'Accuracy CI95 (%)', 'FontSize', font_sz); % overwrite label
colorbar(curr_sp_hdl, 'Location', 'southoutside', ...
    'Ticks', [-1, 0, 1], 'TickLabels', {'< chance', 'null', '> chance'});
set(curr_sp_hdl, 'CLim', [-1.5 1.5]);
title(subplot_titles{3}, 'Interpreter', 'none');

% CI95-ACC relation
sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
% -------------- %
plot_accuracy_CI(result, curr_sp_hdl, do_MCC);
% -------------- %
title('D. Accs.-CI95 relation'); % overwrite title

% Diagonal
sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
% -------------- %
plot_timecourse_CI_filled(result, customticks, customlabels, do_MCC, 1, [], curr_sp_hdl);
% -------------- %
title('E. Diag. time-course of A + Bonf. CI95'); % overwrite title

% Training data slices
sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
% -------------- %
plot_timecourse_CI_filled(result, customticks, customlabels, do_MCC, 0, task_rest_tp, curr_sp_hdl);
% -------------- %
title('F. Training set test accs. + Bonf. CI95'); % overwrite title

%% ROI fig

try
    sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
    c = struct; c.data_type = datatype;
    
    % Check mask file & plot ------------------ %
    [fig_title, maskfile] = check_mask_plotting(result);
    fig_title = plot_maskfile(maskfile, c, fig_title, [], [], curr_sp_hdl);
    % ----------------------------------------- %
    title([sprintf('%s. ', char(64 + sp_ind)) fig_title], 'Interpreter', 'none');
catch e
    e %#ok<NOPRT>
    text(1, 1, '    Retrieving masks failed')
end

%% MCC 2 to 6 (results mostly equal)

for MCC_ind = 2:6
    
    [subplot_stat, subplot_titles, ~, ~] = select_MCC_type(MCC_ind, result, targetfile); % Result/MCC type
    
    if ~isempty(subplot_stat)
        
        if MCC_ind == 2 % all other MCC results use uncorrected CI95, because there are no simple CI methods for permutations or rankings
            
            % CI95
            sp_ind = sp_ind + 1; subplot(n_rows, n_cols, sp_ind); % Sp count
            imagesc(subplot_stat{2}(1:43, 1:43));
            bar_h = decoding_plot_settings(curr_sp_hdl, datatype, font_sz,  customticks, customlabels); % Standard settings
            title(replaceBetween(subplot_titles{2},1,3, sprintf('%s. ', char(64 + sp_ind))), 'Interpreter', 'none'); % replace letter
            xlabel(bar_h, 'Accuracy CI95 (%)', 'FontSize', font_sz);
            
            % CI95-ACC relation
            sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
            [~, plot_title] = plot_accuracy_CI(result, curr_sp_hdl, MCC_ind);
            title(replaceBetween(plot_title,1,3, sprintf('%s. ', char(64 + sp_ind))), 'Interpreter', 'none'); % replace letter
        end
        
        if numel(subplot_stat) == 4 % Multiple more stats eg with TFCE
            
            % P pos
            sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
            imagesc(subplot_stat{2}(1:43, 1:43));
            decoding_plot_settings(curr_sp_hdl, datatype, font_sz,  customticks, customlabels); % Standard settings
            title(replaceBetween(subplot_titles{2},1,3, sprintf('%s. ', char(64 + sp_ind))), 'Interpreter', 'none'); % replace letter
            
            
            % P neg
            sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
            imagesc(subplot_stat{4}(1:43, 1:43));
            decoding_plot_settings(curr_sp_hdl, datatype, font_sz,  customticks, customlabels); % Standard settings
            title(replaceBetween(subplot_titles{4},1,3, sprintf('%s. ', char(64 + sp_ind))), 'Interpreter', 'none'); % replace letter
        end
        
        % H
        sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
        imagesc(subplot_stat{3}(1:43, 1:43));
        bar_h = decoding_plot_settings(curr_sp_hdl, datatype, font_sz,  customticks, customlabels); % Standard settings
        colormap(curr_sp_hdl, [teal; 1 1 1; orange]);
        title(replaceBetween(subplot_titles{3},1,3, sprintf('%s. ', char(64 + sp_ind))), 'Interpreter', 'none'); % replace letter
        xlabel(bar_h, 'Accuracy CI95 (%)', 'FontSize', font_sz);
        colorbar(curr_sp_hdl, 'Location', 'southoutside', ...
            'Ticks', [-1, 0, 1], 'TickLabels', {'< chance', 'null', '> chance'});
        set(curr_sp_hdl, 'CLim', [-1.5 1.5]);
        
        % Diagonal
        sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
        plot_timecourse_CI_filled(result, customticks, customlabels, MCC_ind, 1, [], curr_sp_hdl);
        diag_title = replaceBetween('X. Diag. time-course of A + uncorr. CI95',1,3, sprintf('%s. ', char(64 + sp_ind)));
        title(diag_title, 'Interpreter', 'none'); % replace letter
        
        % Training data slices
        sp_ind = sp_ind + 1; curr_sp_hdl = subplot(n_rows, n_cols, sp_ind); % Sp count
        plot_timecourse_CI_filled(result, customticks, customlabels, MCC_ind, 0, task_rest_tp, curr_sp_hdl);
        slice_title = replaceBetween('X. Training set test accs. + uncorr. CI95',1,3, sprintf('%s. ', char(64 + sp_ind)));
        title(slice_title, 'Interpreter', 'none'); % replace letter

    end
end

%% SAVE

dispv(1, '    %0i statistics plotted for group-level decoding', n_subplots);
save_plot(targetfile, resfig, 0, 0, 0);
% print('-bestfit', resfig, '-dpdf', [targetfile '.pdf']); % as PDF

end