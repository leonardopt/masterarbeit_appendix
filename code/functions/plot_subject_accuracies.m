% function acc_fighdl = plot_subject_accuracies(cfg, ACCs, fig_title, maskfile, layer_descr)
%
% Plots subject-specific accuracies for all comparisons. Also attempts to
% plot a 2D or 3D mesh outline of the mask used for data generation. Needs
% maskfile for that, and layer_descr as well if data is simulated.
%
% Ingmar & Kai, May-2020

% Separated from afx_tdt_subject_timeresolved_decoding.m
% Ingmar 17-03-21

function acc_fighdl = plot_subject_accuracies(cfg, ACCs, fig_title, maskfile, layer_descr)

acc_fighdl = figure('name', fig_title);

for comp_ind = 1:(numel(ACCs) + 1)
    sph = subplot(1, numel(ACCs) + 1, comp_ind);
    
    if (numel(ACCs) + 1) ~= comp_ind
        
        % color img
        imagesc(ACCs{comp_ind}.ACC, [0 100]);
        axis equal; axis tight;
        colormap(customcolormap_preset('pasteljet'));
        cbhdl = colorbar('Location', 'southoutside');
        
        if comp_ind == numel(ACCs)
            ylabel(cbhdl, 'Decoding accuracy (%)', 'Interpreter', 'none', 'FontSize', 10);
        end
        title(ACCs{comp_ind}.title, 'Interpreter', 'none', 'FontSize', 10);
        xlabel('Test time', 'FontSize', 10);
        if comp_ind == 1, ylabel('Train time' , 'FontSize', 10); end
        
        % add tick description
        axislabels = ACCs{comp_ind}.axislabel;
        tick_list = unique(round(linspace(1, length(axislabels.testtime_xaxis_str), 10))); % <--- select a suitable number of ticks
        set(gca, 'XTickLabelRotation', 60, 'FontSize', 10);
        xticks(tick_list); xticklabels(axislabels.testtime_xaxis_str(tick_list));
        yticks(tick_list); yticklabels(axislabels.traintime_yaxis_str(tick_list));
    else
        % plot ROI as last plot
        try
            maskname = ACCs{comp_ind}.maskname;
        catch
            maskname = ''; % avoid problems with real data, maskname only needed for simulation
        end
        if ~exist('layer_descr', 'var')
            layer_descr = '';  % avoid problems with real data, layer_descr only needed for simulation
        end
        fig_title = plot_maskfile(maskfile, cfg, fig_title, maskname, layer_descr);
    end
end % comp

try
    ht = suptitle(fig_title);
    set(ht, 'Interpreter', 'none');
catch
    title(fig_title, 'Interpreter', 'none');
end

drawnow; % reset figure size to show filenames later

end