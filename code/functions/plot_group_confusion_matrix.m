% function result.files = plot_confusion_matrix_results(mean_CONFMAT_els,
% result, result.files, targetfile, plotting_type, customticks, customlabels)
%
% Constructs subplots that show the percentage accurate and missed
% classifications for all label comparisons.
%
% Ingmar, 08-07-20

function result = plot_group_confusion_matrix(mean_CONFMAT_els, result, targetfile, plotting_type, customticks, customlabels)
%% Plotting confusion matrix timeseries for all elements x train x test

figtitle = sprintf('Time-resolved confusion matrix\n%s (mask %0i) with total n = %0i', ...
    result.mask_name, result.mask_index, length(result.subs_todo));
tp_end      = length(result.axislabels.testtime_xaxis_str);
offset      = 0.5; % Axes move 0.5 point automatically because of reasons?
reststart   = 27 + offset;

% Calculate size
conf_mat_sz = sqrt(length(mean_CONFMAT_els)); % Gets number of cols/rows
if mod(conf_mat_sz, 3) ~= 0
    undecided = 1; % Extra undecided column present, otherwise confusion matrix would be square
    if numel(mean_CONFMAT_els) == 12
        nrows = 3; ncols = 4;
    end
else
    nrows = conf_mat_sz; ncols = conf_mat_sz;
end    

% Fontsize
font_sz = 8;
    
%% Standard subplot, done like the subject specific plots

if any(strcmp(plotting_type, 'subplots'))
    cfg.design.conditions2set_map_axislabel = customlabels;
    % -------------------------------- %
    confmat_fighdl = plot_confmats(cfg, mean_CONFMAT_els, undecided, 'All data');
    % -------------------------------- %
    save_plot(confmat_fighdl) % Save
end

%% Concatenate elements into single big matrix

if any(strcmp(plotting_type, 'single_matrix'))
    
    % Initialize figure
    conf_fighdl = figure('name', figtitle);
    targetfile_matrix = [targetfile '_single-matrix'];
    
    % Reshape and plot
    catmat = cat(1, mean_CONFMAT_els{:});     
    orgmat_nrows = size(catmat, 1);
    orgmat_ncols = size(catmat, 2);
    newmat_nrows = orgmat_nrows / ncols;
    newmat_ncols = orgmat_ncols * ncols;
    newmat = nan(newmat_nrows, newmat_ncols);
    for colc_ind = 1:ncols
        newmat_startcol = 1 + (colc_ind-1) * orgmat_ncols;
        newmat_endcol = newmat_startcol + orgmat_ncols - 1;
        catmat_startrow = 1 + (colc_ind-1) * newmat_nrows;
        catmat_endrow = catmat_startrow + newmat_nrows - 1;
        newmat1 = catmat(catmat_startrow:catmat_endrow, :);
        newmat(:, newmat_startcol:newmat_endcol) = newmat1;
    end
    imagesc(newmat);    
    
    % Axes, ticks, colors
    set(gca, 'XTickLabelRotation', 60, 'FontSize', 10);
    [magma, inferno, plasma, viridis] = imported_colormaps; %#ok<*ASGLU>
    colormap(plasma); cbhdl = colorbar; 
    ylabel(cbhdl, 'Classification ratio (%)', 'FontSize', font_sz)
    title(sprintf('%s\n', figtitle), 'Interpreter', 'None', 'FontSize', font_sz);
    set(gca, 'TickDir', 'out'); axis equal; axis tight;
    
    if conf_mat_sz == 3 % 3x3
        
        set(gcf, 'Position', [0 0 800 600]);
        
        % Add separators
        xline(tp_end+offset); xline(90+ offset); yline(tp_end+offset); yline(90+ offset);
        xline(reststart, 'k--'); xline(reststart + tp_end, 'k--'); xline(reststart + tp_end*2, 'k--'); % start of rest phase
        yline(reststart, 'k--'); yline(reststart + tp_end, 'k--'); yline(reststart + tp_end*2, 'k--');
        
        % Labels
        yticks([customticks+ offset (tp_end + offset+ customticks) (tp_end*2 + offset+ customticks)]);
        xticks([customticks+ offset (tp_end + offset+ customticks) (tp_end*2 + offset+ customticks)]);
        yticklabels([customlabels customlabels customlabels]);
        xticklabels([customlabels customlabels customlabels]);
        
        % One way of setting subcategories...
        ylabel(sprintf('Expected labels\nSimilarity      Addition      Congruency'), ...
            'FontSize', font_sz)
        xlabel(sprintf('Congruency      Addition      Similarity\nTrue labels'), ...
            'FontSize', font_sz)
        
    elseif conf_mat_sz == 2 % 2x2
        
        set(gcf, 'Position', [0 0 500 400]);
        
        % Add separators
        xline(tp_end+ offset); yline(tp_end+ offset);
        xline(reststart, 'k--'); xline(reststart + tp_end, 'k--');  % start of rest phase
        yline(reststart, 'k--'); yline(reststart + tp_end, 'k--');
        
        % Labels
        yticks([customticks+ offset (tp_end + offset+ customticks)]);
        xticks([customticks+ offset (tp_end + offset+ customticks)]);
        ylabel(sprintf('Expected labels\nLabel 1        Label 2'), ...
            'FontSize', font_sz)
        xlabel(sprintf('Label 1        Label 2\nTrue labels'), ...
            'FontSize', font_sz)
        
    elseif nrows == 3 && ncols == 4
        
        set(gcf, 'Position', [200 200 850 650]);
        
        % Add separators
        xline(tp_end+offset); yline(tp_end+offset); 
        xline(tp_end*2 +offset); yline(tp_end*2 +offset);
        xline(tp_end*3 +offset);
        xline(reststart, 'k--'); xline(reststart + tp_end, 'k--'); 
        xline(reststart + tp_end*2, 'k--'); xline(reststart + tp_end*3, 'k--');% start of rest phase
        yline(reststart, 'k--'); yline(reststart + tp_end, 'k--'); 
        xline(reststart + tp_end*3, 'k--');
        
        % Labels
        yticks([customticks+ offset (tp_end + offset+ customticks) (tp_end*2 + offset+ customticks)]);
        xticks([customticks+ offset (tp_end + offset+ customticks) (tp_end*2 + offset+ customticks) (tp_end*3 + offset+ customticks)]);
        yticklabels([customlabels customlabels customlabels]);
        xticklabels([customlabels customlabels customlabels customlabels]);
        
        % One way of setting subcategories...
        ylabel(sprintf('True class\nCongruency        Addition        Similarity'), ...
            'FontSize', font_sz)
        xlabel(sprintf('Similarity        Addition        Congruency        Undecided\nPredicted label'), ...
            'FontSize', font_sz)
        
    end
    
    % Save
    save_plot(targetfile_matrix, conf_fighdl, 1, 0, 1)
end

% Custom colormap: heatmap (medium saturation):
% near white - orange yellow - dark pink - dark purple
% confmat_colormap = customcolormap([0 0.25 0.6 1], ...
%    [1 244/255 222/255; 1 0.84 0.21; 237/255 111/255 164/255; 90/255 18/255 150/255]);
