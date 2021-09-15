% function [resfig, targetfile_fig] = plot_decoding_results(result, targetfile, customticks, customlabels, datatype, do_MCC)
%
% Constructs subplots that show hypothesis testing, p-values, mean group
% accuracies and confidence intervals for all training x testing
% timepoints for group condition 1 vs group condition 2 t-test results.
%
% This function only plots. The calculation is done in afx_grouplevel_decoding_stats
%
% Multiple multiple comparison correction options:
%  do_MCC: can be a list with multiple options, e.g. do_MCC = 0:5 for all
%       0: uncorr
%       1: bonf.
%       2: FDR 0.001
%       3: signperm_test clustersize reconstructed from individual quandrant matrices
%       4: signperm_test clustermass reconstructed from individual quandrant matrices
%       5: TFCE
%       6: TFCE reconstructed from individual quandrant matrices
%
% if do_MCC is empty or not existing, a default is set (see function below)
%
% Old custom colormaps for thesis:
% main colormap: '#e9357f', '#ffac1b', '#ffe592', '#b1dee9', '#2b62c0'
% colorstr = {'#2b62c0', '#e9357f', '#ffe592', '#b1dee9', '#2b62c0'};
% colormap(curr_subplot, customcolormap(linspace(0, 1, length(colorstr)), colorstr));
%
% OUT
%   resfig{do_MCC_pos}: handle to figure for test at position do_MCC (not for the do_MCC value)
%   targetfile_fig{do_MCC_pos}: file for resfig{do_MCC_pos}
%
% Kai & Ingmar, May 2020

function [resfig, targetfile_fig] = plot_decoding_results(result, targetfile, customticks, customlabels, datatype, do_MCC)

% Empirical or simulated?
if ~exist('datatype', 'var') || isempty(datatype)
    datatype = 'empirical'; disp('Assuming real data');
end

if ~exist('do_MCC', 'var') || isempty(do_MCC)
    disp('Types of 2nd levels to plot are not specified in do_MCC, using defaults')
    do_MCC = 1:6; % 1:6: everything except uncorr
end

% No group name
if ~isfield(result, 'group_resname')
    result.group_resname = 'Group-level decoding results';
end

%% Plot corrections
% if nothing changed in subplot_creation
% 0: uncorr
% 1: Bonferroni corrected
% 2: FDR
% 3: Sign permutation (clustersize) reconstructed from individual quandrant matrices
% 4: Sign permutation (clustermass) reconstructed from individual quandrant matrices
% 5: TFCE
% 6: TFCE reconstructed from individual quandrant matrices

for mcc_ind = 1:length(do_MCC)
    curr_MCC = do_MCC(mcc_ind);
    [resfig{mcc_ind}, targetfile_fig{mcc_ind}] = subplot_creation(result, customticks, customlabels, curr_MCC, datatype, targetfile); %#ok<AGROW>
    if ~isempty(targetfile_fig{mcc_ind})
        save_as_fig = 1;
        save_plot(targetfile_fig{mcc_ind}, resfig{mcc_ind}, '', save_as_fig);
    else
        fprintf('Cannot save anything for do_MCC = %i because subplot_creation returned an empty resultfig (maybe data was not computed)\n', curr_MCC)
    end
end

end

%% Call subplot creation

function [resfig, targetfile] = subplot_creation(result, customticks, customlabels, curr_MCC, datatype, targetfile)
%% Create figure

% Rows, cols, fig size, font size
n_rows = 3; n_cols = 4; n_subplots = 9; %n_rows * n_cols;
resfig = []; % init empty if return on error
% Plot size
% fig_sz = [300 300 2100 800] % for 1x5 subplots
% fig_sz = [600 0 1315 1218]; % for 2x3 subplots
% fig_sz = [ 30          21        1522         844]; % for 2x4 subplots incl. ROI mask
fig_sz = [63 -353 1169 1243]; % for 3x4 subplots incl. ROI mask
font_sz = 8;

% Which MCC type
[subplot_stat, subplot_titles, MCC_type, targetfile] = select_MCC_type(curr_MCC, result, targetfile);

% Colors
teal = [43 147 166]; teal = teal/255;% lightblue = [53/255 119/255 233/255];
orange = [255 156 62]; orange = orange/255;

% For some reason the axes always offset 0.5 with axes tight
offset = 0.5;

% Init
resfig = figure('name', [result.group_resname, MCC_type], 'Position', fig_sz);

%% Make subplots

sph = gobjects(1, n_subplots); % pre-allocate
for sp_ind = 1:n_subplots
    
    sph(sp_ind) = subplot(n_rows, n_cols, sp_ind);
    if any(sp_ind == 1:4)
        if length(subplot_stat) < sp_ind || isempty(subplot_stat{sp_ind})
            % skip step if not data is present
            try delete(sph(sp_ind)); end % remove initiated subplot,
            continue
        end
        
        shading('flat');
        curr_stat = squeeze(subplot_stat{sp_ind});
        imagesc(sph(sp_ind), curr_stat(1:43, 1:43));
        
        if any(contains(subplot_titles(sp_ind), ...
                {'significance', 'hypothesis', 'test'}, 'IgnoreCase', true))
            
            % colorbar not needed for hypothesis testing
            %             if any(subplot_stat{sp_ind} == -1, 'all') % below-chance accs
            colormap(sph(sp_ind), [teal; 1 1 1; orange]);
            colorbar(sph(sp_ind), 'Location', 'southoutside', ...
                'Ticks', [-1, 0, 1], 'TickLabels', {'< chance', 'null', '> chance'});
            %             elseif ~any(subplot_stat{sp_ind} == 1, 'all')
            %                 colormap(sph(sp_ind), [1 1 1]); % no effect at all
            %                 colorbar(sph(sp_ind), 'Location', 'southoutside', 'Ticks', 0, 'TickLabels', {'null'});
            %             else
            %                 colormap(sph(sp_ind), [1 1 1; orange]); % only above-chance
            %                 colorbar(sph(sp_ind), 'Location', 'southoutside', ...
            %                 'Ticks', [0, 1], 'TickLabels', {'null', '> chance'});
            %             end
            
            % set clims for hypothesis testing -1 to 1
            set(sph(sp_ind), 'CLim', [-1.5 1.5]);
            
            cline = [0 0 0]; % [0.6 0.6 0.6];
            xline(27 + offset, 'Color', cline, 'LineWidth', 1); % 'LineStyle', '--');
            yline(27 + offset, 'Color', cline, 'LineWidth', 1); % 'LineStyle', '--');
        else
            
            bar_h = decoding_plot_settings(datatype, font_sz,  customticks, customlabels);
            
            % Alter color limits
            if isfield(result, 'extend_caxis') && ~isempty(result.extend_caxis)
                if contains(subplot_titles(sp_ind), 'width')
                    caxis(result.extend_caxis(2, :))
                else
                    caxis(result.extend_caxis(1, :))
                end
            end
            
            % add xlabel for accuracies
            if contains(subplot_titles(sp_ind), 'ccurac')
                xlabel(bar_h, 'Decoding accuracy (%)');
            else
                xlabel(bar_h, 'Accuracy CI95 (%)', 'FontSize', font_sz);
            end
        end
        
        % Get axes ticks and labels
        title(subplot_titles{sp_ind}, 'Interpreter', 'none');
        xlabel('Testing time-points'); ylabel('Training time-points');
        yticks(customticks + offset); xticks(customticks + offset);
        yticklabels(customlabels); xticklabels(customlabels);
        
        if strfind(subplot_titles{sp_ind}, 'P_neg')
            set(gca, 'clim', [0 1]);
        end
        
    elseif sp_ind == 4+1 % Diagonal timeslice
        % -------------- %
        sph(sp_ind) = plot_timecourse_CI_filled(result, customticks, customlabels, ...
            curr_MCC, 1, [], sph(sp_ind));
        % -------------- %
        
    elseif  sp_ind == 4+2 % Training result slices
        task_rest_tp = [14 14];
        % -------------- %
        sph(sp_ind) = plot_timecourse_CI_filled(result, customticks, customlabels, ...
            curr_MCC, 0, task_rest_tp, sph(sp_ind));
        % -------------- %
        
    elseif sp_ind == 4+3  % CI95-ACC result relation
        % -------------- %
        sph(sp_ind) = plot_accuracy_CI(result, sph(sp_ind), curr_MCC);
        % -------------- %
        
    elseif sp_ind == 4+4 % ROI
        disp('    Trying to plot ROI, should work for empirical data, not sure about simulation')
        try
            % Check mask file -------- %
            [fig_title, maskfile] = check_mask_plotting(result);
            % ------------------------ %
            
            % Plot mask results ------ %
            c.data_type = datatype;
            fig_title = plot_maskfile(maskfile, c, fig_title);
            xlabel(fig_title, 'Interpreter', 'none');
            % ------------------------ %
        catch e
            e %#ok<NOPRT>
            text(0, 0, 'retrieving masks failed')
        end
        
    elseif sp_ind == 4 + 5 % Accuracy contour for all possible MCC options above
        
        sph(sp_ind) = subplot(n_rows, n_cols, sp_ind:sp_ind+1); % enlarge
        % -------------- %
        sph(sp_ind) = plot_MCC_contours(result, sph(sp_ind), customticks, customlabels);
        % -------------- %
        
    else
        warning('Unknown subplot sp_ind %i, perhaps empty space.', sp_ind)
    end
    try % try to make text + title smaller
        %         h = get(gca, 'Title');
        %         h.FontSize = 8;
        set(gca, 'FontSize', 8);
        set(gca, 'TitleFontSizeMultiplier', 1.1);
    catch e
        e %#ok<NOPRT>
        disp('Could not change FontSize or TitleFontSizeMultiplier in subplot')
    end
end

% add a title for the full figure (take it's name)
try
    h = suptitle(get(gcf, 'Name'));
    set(h, 'Interpreter', 'none');
    set(h, 'FontWeight', 'bold');
end

dispv(1, '    %0i statistics plotted for group-level decoding', n_subplots);
end