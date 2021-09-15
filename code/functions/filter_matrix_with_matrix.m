% Function to filter values in a matrix using another matrix as filter 
% Optional argument 'plot' to plot result 


function filtered_mat = filter_matrix_with_matrix(matrix_to_filter, mat_filter, plot)

mat_filter = logical(mat_filter); % make sure it's a logical array 

% Initialise filtered matrix 
filtered_mat = nan(size(matrix_to_filter));

% Populate filtered matrix with values from the matrix to filter, using the
% filter 
filtered_mat(mat_filter) = matrix_to_filter(mat_filter);

if exist('plot', 'var')
    
    masknans = ones(size(filtered_mat));
    masknans(isnan(filtered_mat)) = 0;
    imagesc(filtered_mat, 'AlphaData', masknans)
    colorbar
    title('Selected values (white where NaN)')
    
end

end