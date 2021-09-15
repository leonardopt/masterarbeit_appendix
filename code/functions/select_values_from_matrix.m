function mat = select_values_from_matrix(filterrows, filtercolumns, plot)

% Make sure filterrows is a column vector and filtercolumns is a row vector 
% filterrows will filter the rows in the matrix, filtercoulmns the columns

if size(filterrows, 2) ~= 1
    filterrows = filterrows';
    
end

if size(filtercolumns, 1) ~= 1
    filtercolumns = filtercolumns';
end 


mat = and(repmat(filterrows,[1,size(filterrows, 1)]), repmat(filtercolumns,[size(filtercolumns, 2), 1]));

if exist('plot', 'var')
figure
imagesc(mat)
title('Selected values')
ylabel(inputname(1))
xlabel(inputname(2))
end 

end