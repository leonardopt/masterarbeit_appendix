% Function that creates a filter to get a target string in a cell array

function filter = filter_condition(source_labels, target_string1, target_string2)

filter = contains(source_labels, target_string1);
if exist('target_string2', 'var') 
    filter2 = contains(source_labels, target_string2);
    filter = and(filter, filter2);     
end 

end