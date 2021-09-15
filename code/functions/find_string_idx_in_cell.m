function position_idx = find_string_idx_in_cell(all_cells, string_to_find)


position_idx = find(~cellfun(@isempty, regexp(all_cells, string_to_find)));


end 