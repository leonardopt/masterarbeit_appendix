% function [newmat, newlabelnames, newlabelnames2] = average_mat(mat, index, labelnames, index2, labelnames2, mean_fun)
%
% This functions builds average matrices for combinations of rows and
% colums that are passed as "index" and "index2". See example below how it 
% works. To know which rows and columns are averaged in the output matrix,
% a stringcell input (labelnames, labelnames2) with row and column names is 
% used to create new labelnames for the output matrix.
%
% If mat has a 3rd dimension, the newmat will also have 3 dimensions and
% averaging is done for the first 2 (these are assumed images of
% indidividual subjects).
%
% mean_fun (default @mean) sets the function to use to calculate the mean
% (or any other summary you wish). Use e.g. @nanmean if nans are contained 
% in the data that you like to ignore.
% 
% 
%
% Example code
% %% Example data
% mat = [11 12 13
%     21 22 23
%     31 32 33]
% labelnames = {'r1', 'r2', 'r3'} % rows
% labelnames2 = {'c1', 'c2', 'c3'} % columns
% 
% %% Example calls
% % average all data, i.e. put all rows and colums in new class "1"
% index = [1 1 1]; % values are not important except that the result will be ordered increasing
% % if only 1 index is passed it will be used for rows and columns
% [newmat, newlabelnames] = average_mat(mat, index, labelnames);
% 
% % average for each column all rows
% index = [4 4 4]; % for rows, Again, exact values are not important.
% index2 = [4 5 6]; % for columns. values of index and index2 are not related in any way
% [newmat, newlabelnames, newlabelnames2] = average_mat(mat, index, labelnames, index2, labelnames2);
% 
% % average upper right(1:2, 2:3) and lower left (1, 3)
% index = [1 1 2]; % for rows, Again, exact values are not important.
% index2 = [2 1 1]; % for columns. values of index and index2 are not related in any way
% [newmat, newlabelnames, newlabelnames2] = average_mat(mat, index, labelnames, index2, labelnames2);
% 
% % ignore row and colum 2 and average the rest
% % average all data, i.e. put all rows and colums in new class "1"
% index = [1 nan 1]; % values are not important except that the result will be ordered increasing
% % if only 1 index is passed it will be used for rows and columns
% [newmat, newlabelnames] = average_mat(mat, index, labelnames);
%
% % select lower left without averaging (the same as mat(2:3, 1:2)
% index = [1 2 nan]; % for rows, Again, exact values are not important.
% index2 = [nan 1 2]; % for columns. values of index and index2 are not related in any way
% [newmat, newlabelnames, newlabelnames2] = average_mat(mat, index, labelnames, index2, labelnames2);
% Kai, 2017/01/12 (adapted from v2017/01/12)

% History: 
% - sorting of indices removed 
% Leonardo 09.08.2019
%

% Potential improvement: nan handling in index vectors
% Taking out rows and columns using nan is
% implemented very ineffictive and confusing, could be done much simpler
% e.g. by first detecting if there are nans as index 1 or 2, then running
% the function again without these rows/columns, and then returning the
% result


function [newmat, newlabelnames, newlabelnames2] = average_mat(mat, index, labelnames, index2, labelnames2, mean_fun)

% Defaults
if ~exist('index2', 'var') || isempty(index2)
    index2 = index;
end
if ~exist('mean_fun', 'var') || isempty(mean_fun)
    mean_fun = @mean;
end

% Get unique elements, keep order
[u_ind,ia] = unique(index);
%u_ind = index(sort(ia)); % sort by first occurence

[u_ind2,ia2]  = unique(index2);
%u_ind2 = index2(sort(ia2)); % sort by first occurence

% store rows and columns to take out
row_out = [];
col_out = [];

% Create mean and new labelnames
for ind = 1:length(u_ind)
    % Continue wit next unique index if current u_ind is isnan
    if isnan(u_ind(ind)), row_out(end+1) = ind; continue, end
    % find all elements that are in the current class
    ind_filt = find(index == u_ind(ind));
    
    if exist('labelnames', 'var')
        newlabelnames{ind} = [labelnames{ind_filt}];
    end
    
    for ind2 = 1:length(u_ind2)
        % Continue wit next unique index if current u_ind is isnan
        if isnan(u_ind2(ind2)), col_out(end+1) = ind2; continue, end
        % find all elements that are in the current class
        ind2_filt = find(index2 == u_ind2(ind2));
        if exist('labelnames2', 'var') && ~isempty(labelnames2)
            newlabelnames2{ind2} = [labelnames2{ind2_filt}];
        end
        for s3_ind = 1:size(mat, 3)
            curr_mat = mat(ind_filt, ind2_filt, s3_ind);
            newmat(ind, ind2, s3_ind) = mean_fun(curr_mat(:)); % default: mean()
        end
    end
end

% Remove all rows and columns that contained nans
% nans at the end have not been added, so we do not need to remove them
row_out = unique(row_out);
row_out(row_out > size(newmat, 1)) = [];
col_out = unique(col_out);
col_out(col_out > size(newmat, 2)) = [];
newmat(row_out, :, :) = [];
newmat(:, col_out, :) = [];
try % also remove label names, in case they exist
    newlabelnames1(row_out) = [];
    newlabelnames2(col_out) = [];
end

