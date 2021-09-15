function mark_nans(m, marker, markercolor, warntext_sprintf, warntextcolor)
% add nans
nan_pos = find(isnan(m));

% defaults

if ~exist('marker', 'var'), marker = 'x'; end
if ~exist('markercolor', 'var'), markercolor = 'r'; end
if ~exist('warntext_sprintf', 'var'), warntext_sprintf = 'Found %i nans in matrix to plot'; end
if ~exist('warntextcolor', 'var'), warntextcolor = 'r'; end

if ~isempty(nan_pos)
    n_nans = length(nan_pos);
    warntext = sprintf(warntext_sprintf, n_nans);
    warning([warntext ', marking it in plot'])
    
    [X, Y] = meshgrid(1:size(m, 2), 1:size(m, 1));
    nan_x = X(nan_pos);
    nan_y = Y(nan_pos);
    hold all
    scatter(nan_x', nan_y', marker, markercolor);
    
    text(1, 1, warntext, 'color', warntextcolor);
    
end

end