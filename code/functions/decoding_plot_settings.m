% bar_h = decoding_plot_settings(ax_hdl, datatype, font_sz,  customticks, customlabels)
%
% Some plot settings for plot_decoding_results.m and similar functions.
% Only result struct is required (see grouplevels and 
% afx_grouplevel_decoding_stats), other arguments are optional.
%
% Ingmar, 10-09-21

function bar_h = decoding_plot_settings(ax_hdl, datatype, font_sz,  customticks, customlabels)

if ~exist('ax_hdl', 'var') || isempty(ax_hdl)
    disp('    decoding_plot_settings: ax_hdl missing as first argument, assuming current axes');
    ax_hdl = gca;
    % fig_hdl = gcf;
end

% Colormap
shading('flat');
[~, ~, plasma, viridis] = imported_colormaps;
if exist('datatype', 'var') && contains(datatype, 'simulat')
    colormap(ax_hdl, plasma);
else
    colormap(ax_hdl, viridis);
end

% Font size
if ~exist('font_sz', 'var')
    font_sz = 8;
end
set(gca,'FontSize', font_sz);
bar_h = colorbar(ax_hdl, 'Location', 'southoutside', 'FontSize', font_sz);

% Round to 2 dec, correct ends of colorbar
bar_h.Ticks   = round(linspace(bar_h.Limits(1), bar_h.Limits(2), 5), 2);
bar_h.Limits  = round(bar_h.Limits, 2);

% Add line
cline = [0 0 0]; % [0.6 0.6 0.6];
offset = 0.5;
xline(27 + offset, 'Color', cline,  'LineWidth', 1);
yline(27 + offset, 'Color', cline,  'LineWidth', 1);

% Visualization
axis equal; axis tight;
set(ax_hdl, 'TickDir', 'out');
set(ax_hdl, 'XTickLabelRotation', 60); 

% Labels and ticks
if ~exist('customticks', 'var') || ~exist('customlabels', 'var')
    [customticks, customlabels] = afx_set_timings(0);
end    
yticks(customticks + offset); xticks(customticks + offset);
yticklabels(customlabels); xticklabels(customlabels);
xlabel('Testing time-points'); ylabel('Training time-points');
xlabel(bar_h, 'Decoding accuracy (%)');

disp('    Subplot defaults set for AFX decoding: decoding_plot_settings.m');
end