% function save_plot(plot_fname, fighdl, info, save_as_fig, save_as_png)
%
% Saves current plot or specified plot handle. On default, saves in .png
% and .svg, and possibly .fig. Separated from subject-specific decoding as
% it could be useful elsewhere. Does not neccessarily need TDT but does try
% to use it first (as wrapper for save_fig.m, not to be confused with 
% savefig from MATLAB).
%
% IN
%    plot_fname: should be a full path including file name, but without
%    file extension.
%
% OPTIONAL:
%    fighdl: figure handle if specified, otherwise current open plot
%    info: if set to 1, try to add date, time and file location to the plot xlabel
%    save_as_fig: if set to 1, will save as MATLAB figure .fig
%    save_as_png: if set to 1, will save as .png; now disabled since nearly
%     all plots require some vector graphics editing and to save space.
%
% Ingmar, June 2020

function save_plot(plot_fname, fighdl, plot_info, save_as_fig, save_as_png)
%% Saving options

% Check if path is valid
try
    [path, name, ext] = fileparts(plot_fname);
    isfolder(path);
    if ~isempty(ext)
        disp('    "plot_fname" was passed with file extension, ignoring it...')
        plot_fname = fullfile(path, name);
    end
catch e1
    warning('Specified filepath seems invalid...');
    disp(['    Error: ' e1.message]);
end

% Add time and filename start & endtime, if available
if exist('plot_info', 'var') && ~isempty(plot_info) && (isscalar(plot_info) || ...
        islogical(plot_info)) && plot_info
    try
        disp(['    Appending datestring & file location to ' name ' plot']);
        if length(path) >= 80, path = [path(1:80), '...']; end
        outtext = [datestr(now), newline path];
        
        % Append to handle if possible
        if exist('fighdl', 'var') && ~isempty(fighdl)
            if isfield(fighdl, 'CurrentAxes') && isfield(fighdl.CurrentAxes, 'Xlabel')
                new_xlabel = get(fighdl.CurrentAxes.XLabel, 'String');
            else
                axeshdl = gca;
                new_xlabel = axeshdl.XLabel.String;
            end
            new_xlabel = [new_xlabel newline newline outtext]; % Append datestr
            xlabel(new_xlabel, 'Interpreter', 'none');
        elseif ~exist('fighdl', 'var')
            new_xlabel = gca('XLabel');
            new_xlabel = [new_xlabel newline newline outtext]; % Append datestr
            xlabel(new_xlabel, 'Interpreter', 'none');
        else
            xlabel([newline outtext], 'Interpreter', 'none', 'FontSize', 6);
        end
    catch e2
        warning('Could not add time and/or filename to figure')
        disp(['    Error: ' e2.message]);
        return;
    end
end

%% Saving

% Change here for other vector graphic formats (eg PDF/EPS)
default_format  = '-dsvg'; 
default_ext     = '.svg';

% fig
if exist('save_as_fig', 'var') && ~isempty(save_as_fig) && ...
        (isscalar(save_as_fig) || islogical(save_as_fig)) && save_as_fig
    plot_cfg.plot_design_formats = {default_format};
else
    plot_cfg.plot_design_formats = {default_format, 'NOFIGFILE'};
end

% png
if exist('save_as_png', 'var') && ~isempty(save_as_png) && ...
        (isscalar(save_as_png) || islogical(save_as_png)) && save_as_png
    plot_cfg.plot_design_formats{end+1} = '-dpng';
end

plot_fname_backup = plot_fname; % For potential warning later
try
    try % With TDT
        if exist('fighdl', 'var')
            save_fig(plot_fname, plot_cfg, fighdl);
        else % Currently open figure
            save_fig(plot_fname, plot_cfg);
        end
        disp(['    Plot saved: ' plot_fname]); layout_line_break(1);
    catch % If TDT is missing
        disp('    save_fig() failed, trying standard MATLAB saving...');
        if exist('fighdl', 'var')
            
            % Default
            plot_fname = [plot_fname, default_ext]; % For potential warning below
            saveas(fighdl, plot_fname)
            
            % Standard PNG
            if exist('save_as_png', 'var') && ~isempty(save_as_png) && ...
                    (isscalar(save_as_png) || islogical(save_as_png)) && save_as_png
                saveas(fighdl, [plot_fname_backup '.png'])
            end
            
            % Optional fig
            if exist('save_as_fig', 'var') && ~isempty(save_as_fig) && ...
                    (isscalar(save_as_fig) || islogical(save_as_fig)) && save_as_fig
                savefig(fighdl, [plot_fname_backup '.fig'])
            end
        else
            % Default
            plot_fname = [plot_fname, default_ext]; % For potential warning below
            saveas(plot_fname)
            
            % Standard PNG
            if exist('save_as_png', 'var') && ~isempty(save_as_png) && ...
                    (isscalar(save_as_png) || islogical(save_as_png)) && save_as_png
                saveas(fighdl, [plot_fname_backup '.png'])
            end
            
            % Optional fig
            if exist('save_as_fig', 'var') && ~isempty(save_as_fig) && ...
                    (isscalar(save_as_fig) || islogical(save_as_fig)) && save_as_fig
                savefig(fighdl, [plot_fname_backup '.fig'])
            end
        end
    end
catch e3
    warning('Plot could not be saved at %s.', plot_fname);
    disp(['Error: ' e3.message]);
end
fprintf('\n    Plot saved: %s\n', plot_fname);  
end