% Function to make labels for the matrices (mean, CI in c/s separate and
% merged conditions)

function [row_labelnames, column_labelnames] = set_labels_for_matrices(row_labelnames, column_labelnames)

set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(0, 'DefaultAxesTickLabelInterpreter', 'none')

%title('Time-resolved correlation between odd and even runs, averaged across runs and participants');
row_labelnames = cellfun(@(x) regexp(x, '.+?(?=[*])', 'once', 'match'), row_labelnames, 'UniformOutput', false);
column_labelnames = cellfun(@(x) regexp(x, '.+?(?=[*])', 'once', 'match'), column_labelnames, 'UniformOutput', false);
t_start_row = find_string_idx_in_cell(row_labelnames, 'Cue  00');
t_start_column = find_string_idx_in_cell(column_labelnames, 'Cue  00');
t_rest_row = find_string_idx_in_cell(row_labelnames, 'Rest 00');
t_rest2_row = find_string_idx_in_cell(row_labelnames, 'Rest 18');
t_rest_column = find_string_idx_in_cell(column_labelnames, 'Rest 00');
t_rest2_column = find_string_idx_in_cell(column_labelnames, 'Rest 18');
t_trial_row = find_string_idx_in_cell(row_labelnames, 'Trial28');
t_trial_column = find_string_idx_in_cell(column_labelnames, 'Trial28');



% If merged complex/simple
if length(row_labelnames)< 200
    
    row_labelnames    = cellfun(@(x) strrep(x, 'CR', '_R'), row_labelnames, 'UniformOutput', false);
    column_labelnames = cellfun(@(x) strrep(x, 'CR', '_R'), column_labelnames, 'UniformOutput', false);
    row_labelnames    = cellfun(@(x) strrep(x, 'C.', '_'), row_labelnames, 'UniformOutput', false);
    column_labelnames = cellfun(@(x) strrep(x, 'C.', '_'), column_labelnames, 'UniformOutput', false);
        


end 


% gap = 18;
% ticks_odd_av = unique([1:gap:length(odd_labelnames), t_start_odd, t_rest_odd, t_rest2_odd]);
% ticks_even_av = unique([1:gap:length(even_labelnames), t_start_even, t_rest_even, t_rest2_even]);
ticks_row_av = unique([ t_start_row, t_rest_row, t_rest2_row,t_trial_row ]);
ticks_column_av = unique([ t_start_column, t_rest_column, t_rest2_column, t_trial_column]);
set(gca, 'YTick', ticks_row_av, 'YTickLabel', row_labelnames(ticks_row_av));
set(gca, 'Xtick', ticks_column_av, 'XTickLabel',  column_labelnames(ticks_column_av),'XTickLabelRotation', 90);

line([0,length(row_labelnames)], [t_start_row(:)- 0.5 ,t_start_row(:)- 0.5 ], 'Color', 'w');
line([t_start_column(:)- 0.5 ,t_start_column(:) - 0.5 ], [0,length(column_labelnames)], 'Color', 'w');

xlabel('Even runs')
ylabel('Odd runs')

end 