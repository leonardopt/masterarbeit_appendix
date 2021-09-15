% Function to select congruent conditions (Aeh-Aeh, Kon-Kon, Sum-Sum)
% either for Trial or Rest (to be specified as third argument)


function filtered_mat = select_incongruent_conditions(source_data, source_labels,  trial_rest, trial_rest2)

conds = {'Aeh', 'Kon', 'Sum'};

filtered_mat = nan(length(source_labels));

if exist('trial_rest2', 'var')
       for c = 1:length(conds)
        curr_cond = conds{c};
        frow  =  ~filter_condition(source_labels, curr_cond) ;
        fcol  =  filter_condition(source_labels, curr_cond);
        fmat = select_values_from_matrix(frow, fcol); % first filter: rows, second filter: columns
        current_filtered_mat = filter_matrix_with_matrix(source_data, fmat);
        % Store in the mat
        filtered_mat(fmat) = current_filtered_mat(fmat);
    end 
    
    
else
    for c = 1:length(conds)
        curr_cond = conds{c};
        frow  =  and(~filter_condition(source_labels, curr_cond) , filter_condition(source_labels, trial_rest));
        fcol  =  filter_condition(source_labels, curr_cond);
        fmat = select_values_from_matrix(frow, fcol); % first filter: rows, second filter: columns
        current_filtered_mat = filter_matrix_with_matrix(source_data, fmat);
        % Store in the mat
        filtered_mat(fmat) = current_filtered_mat(fmat);
    end
end