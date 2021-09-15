% [fighdl, plot_title] = plot_accuracy_CI(result, fighdl, do_MCC)
%
% Visualized how much the accuracy and CI95 distributions overlap, showing
% their (potentially) linear relationship. Merged to be a subplot of
% plot_decoding_results.m
%
% Ingmar, 13-10-20

% Updated as subplot for decoding results - Ingmar 23-07-21

% Colors etc. taken from Leo's afx_tool_plot_mean_vs_CI95
% colors = [ [102, 194, 165]; 
%             [141, 160, 203]; 
%             [252, 141, 98]];

function [fighdl, plot_title] = plot_accuracy_CI(result, fighdl, do_MCC)
%% Prepare

% Which multiple comparison correction (MCC) for CI95
sp_stat = {};
if ~exist('do_MCC', 'var') || isempty(do_MCC) || any(do_MCC == 2:6)
    sp_stat{1} = result.uncorr.CItwos_diff; % uncorrected, FDR, sign perm
elseif exist('do_MCC', 'var') && do_MCC == 1
    sp_stat{1} = result.bonf.CItwosided_diff;
end
sp_stat{2} = result.mean;

plot_title = 'F. Accuracy-CI95 relation';
if ~exist('fighdl', 'var') || isempty(fighdl)
    fighdl = figure('name', plot_title, 'Position', [300 300 600 520]);
end

%% Scatterplot groups

% Specify timing
task_end    = 27;
rest_begin  = 28;
tp_end      = size(sp_stat{1}, 1);

colors = [[53 119 233];
        [72 23 102]%[105 55 168];
        [243 80 122]];
colors = colors/255; % RGB

% Pre-allocate group indices
g_inds = zeros(size(sp_stat{1}));

% Task-task  & task-rest
for s_ind = 1:task_end
    g_inds(s_ind, 1:task_end) = 1;
    g_inds(s_ind, rest_begin:tp_end) = 2; 
end

% Rest-rest % rest-task
for s_ind = rest_begin:tp_end
    g_inds(s_ind, 1:task_end) = 2;
    g_inds(s_ind, rest_begin:tp_end) = 3; 
end

%% Plot

hold on;
for s_ind = unique(g_inds)'
    scatter(sp_stat{1}(g_inds == s_ind), sp_stat{2}(g_inds == s_ind), 2, colors(s_ind, :), '.');
end
hold off;

% Chance line
yline(result.chancelevel, 'k-', 'Chance level');

% Labels and titles
plotnames = {'Trial-trial', 'Trial-rest', 'Rest-rest'};
legend(plotnames, 'Location', 'southoutside');
set(gca, 'FontSize', 8); axis square; 
ylabel('Decoding accuracy (%)'); xlabel('CI95-width (%)');
title(plot_title);

end

% Old additional code

% Redundancy
% if ~isfield(result, 'corrcoef_mat') || ~isfield(result, 'corrcoef_p')
%     [result.corrcoef_mat, result.corrcoef_p, result.corrcoef_low, result.corrcoef_up] = corrcoef(result.mean, ...
%         squeeze(diff(result.bonf.CItwosided)),'Alpha',result.bonf.alpha);
% end
% if numel(result.corrcoef_mat) ~= 1
%     result.corrcoef_mat = result.corrcoef_mat(2);
%     result.corrcoef_p   = result.corrcoef_p(2);
% end

% gscatter seems bugged as it only shows 2 groups instead of 3
%gscatter(subplot_stat{1}, subplot_stat{2}, g_indices, colors, '...', 8, 'off')


% Mean centers
% for s_ind = unique(g_inds)'
%     plot(mean(sp_stat{1}(g_inds == s_ind)), mean(sp_stat{2}(g_inds == s_ind)), 'd');
% end