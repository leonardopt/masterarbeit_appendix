% confmat_fighdl = plot_confmats(cfg, CONFMAT, undecided, name)
%
% Plots subject-specific confusion matrices given all labels.
%
% Ingmar, 2020

function confmat_fighdl = plot_confmats(cfg, CONFMAT, undecided, name)

if ~exist('name', 'var')
    name = sprintf('subject %02i', cfg.subj);
end
figtitle = sprintf('Train-by-test time-point confusion matrix: %s', name);
axislabels  = cfg.design.conditions2set_map_axislabel{:};
try
    tick_list = unique(round(linspace(1, length(axislabels.testtime_xaxis_str), 8))); % <--- select a suitable number of ticks
catch
    [tick_list, axislabels] = afx_set_timings(1);
end
if ~undecided
    CONFMAT_plot = cat(1, CONFMAT{:}); % "flatten" cell array to struct
    confmat_fighdl = figure('name', figtitle, 'Position', [300  300  800  800]);
else
    confmat_fighdl = figure('name', figtitle, 'Position', [300  300  1050  800]);
    try
        CONFMAT_plot = [CONFMAT(1,1:4), CONFMAT(2,1:4), CONFMAT(3,1:4)]; % if unflattened
    catch
        CONFMAT_plot = CONFMAT;
    end
end
rowc = 3; colc = 3;

for sp_ind = 1:length(CONFMAT_plot)
    if undecided
        colc = 4;
        curr_subplot = subplot(rowc, colc, sp_ind, 'align');
        if ~isfield(CONFMAT_plot{1}, 'MAT')
            CONFMAT_to_plot{sp_ind}.MAT = CONFMAT_plot{sp_ind};
            clear CONFMAT_plot;
            CONFMAT_plot{sp_ind}.MAT = CONFMAT_to_plot{sp_ind}.MAT;
        end
        imagesc(curr_subplot, CONFMAT_plot{sp_ind}.MAT, [0 100])
        try
            title(curr_subplot, CONFMAT_plot{sp_ind}.compstr, 'Interpreter', 'none', 'FontSize', 6);
        catch
            title(curr_subplot, sprintf('Element %i', sp_ind), 'Interpreter', 'none', 'FontSize', 6);
        end
    else
        curr_subplot = subplot(rowc, colc, sp_ind, 'align');
        imagesc(curr_subplot, CONFMAT_plot(sp_ind).MAT, [0 100])
        title(curr_subplot, CONFMAT_plot(sp_ind).compstr, 'Interpreter', 'none', 'FontSize', 6);
    end
    set(curr_subplot, 'XTickLabelRotation', 60);
    
    % Axes, ticks, colors
    [~, ~, plasma, ~] = imported_colormaps;
    colormap(plasma);
    axis equal; axis tight; set(curr_subplot, 'TickDir', 'out');
     
    % Only yticks on leftmost side & xticks on bottom row
    if any(sp_ind == 1:colc:length(CONFMAT_plot))
        yticks(curr_subplot, tick_list);
        yticklabels(curr_subplot, axislabels.traintime_yaxis_str(tick_list));
        ylabel('Train');
    else
        set(curr_subplot, 'yticklabel', []);
    end
    if any(sp_ind == (colc + colc + 1):length(CONFMAT_plot))
        xticks(curr_subplot, tick_list);
        xticklabels(curr_subplot, axislabels.testtime_xaxis_str(tick_list));
        xlabel('Test');
    else
        set(curr_subplot, 'xticklabel', []);
    end
end

try % Adjust fig position and add 1 colorbar for all subplots
    cbpos = get(subplot(rowc, colc, sp_ind), 'Position');
    chdl = colorbar('Position', [cbpos(1)+cbpos(3) + 0.07 cbpos(2) 0.1 cbpos(2)+cbpos(3) * 3], ...
        'AxisLocation','in');
    title(chdl, '%');
catch
    dispv(1, '    Adding colorbar to CONFMAT display failed...');
end

% Add title and path to fig
suptitle(figtitle); drawnow;

end