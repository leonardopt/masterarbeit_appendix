% [fig_title, phdl] = plot_maskfile(maskfile, cfg, fig_title, maskname, layer_descr, sph)
%
% Plots a maskfile
%
% IN
%   maskfile: currently filename, could be changed to also take mask as matrix directly
%   cfg: struct with fields
%       cfg.data_type: if not string 'simulated', will be plotted as is
%                      if is string 'simulated', try to add extra info
% OPTIONAL
%   fig_title: title for figure (could be made optional)
%   for 'simulated' only
%   maskname: name of mask as string, e.g. from 
%   layer_descr: description of layers
%   sph: subplot handle
%
% OUT
%   title possibly appended with sim mask layer
%   phdl: mesh render object of mask
%
% Ingmar, 2020

% Extracted from plot_subject_accuracies, Kai 2021-09-07
% Kept aspect ratios of axes for consistent mask renders, Ingmar 9-9-21

function [fig_title, phdl] = plot_maskfile(maskfile, cfg, fig_title, maskname, layer_descr, sph)
if ~exist('fig_title', 'var')
    fig_title = 'Mask';
end

% get mask coverage
curr_mask_mat = spm_read_vols(spm_vol(maskfile));

% Add figure with mask coverage
if ~strcmp(cfg.data_type, 'simulated') % 3D brain
    phdl = patch(isosurface(curr_mask_mat, 'noshare')); % plots mesh
    phdl.FaceColor = 'red'; phdl.EdgeColor = 'none';
    set(gca, 'XLim', [0 64], 'YLim', [0 64], 'ZLim', [0 33]);
    pdhl.DataAspectRatioMode = 'manual'; %#ok<STRNU> % seems to keep aspect ratio
    camup([0 0 1]); campos([400 -155 165]); 
    fig_title = [fig_title, ' 3D render'];
    title(fig_title);
else % sim voxels
    % get index of
    [~, ~, dim3] = ind2sub(size(curr_mask_mat), find(curr_mask_mat == 1));
    imagesc(squeeze(curr_mask_mat(:, :, dim3(1)))); % plots 2d img
    axis equal; axis tight; yticks(1:5); xticks(1:5);
    colorbar('Ticks',[0.5, 1.5], 'TickLabels', {'Mask-out', 'Mask-in'}, 'Location', 'southoutside');
    colormap(sph, [1 1 1; 0 0 0]);

    try  % add info in title
        if contains(maskfile, maskname, 'IgnoreCase', true) % check mask name

            % assumes layer number at last pos
            layer = str2double(extractAfter(maskname, strlength(maskname) - 1));
            if ~isempty(layer_descr) && layer == dim3(1)
                fig_title = [fig_title '_' layer_descr(layer)]; 
            end
        end
    catch
        dispv(1, '    Possibly, not all information on sim. mask and/or parameter set was added to plot.');
    end
    
end