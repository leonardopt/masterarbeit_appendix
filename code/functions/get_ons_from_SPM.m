% function [u, ons, conds, filters] = get_ons_from_SPM(SPM)
%
% Function that extracts condition names and onsets by condition from an
% SPM.mat. Moreover, it returns the filters created to extract each
% condition.

function [u, ons, conds, filters] = get_ons_from_SPM(SPM)

conds   = {'left', 'right', 'aeh-s', 'aeh-c', 'kon-s', 'kon-c', 'sum-s', 'sum-c'};
filters = struct();

% Get condition names by run
for run = 1:length(SPM.Sess)
    
    for cnd = 1:length(conds)
        
        % Get current condition
        curr_cond = conds{cnd};
        curr_cond_field = strrep(curr_cond, '-', '_');
        
        % Extract names of the condition from the SPM mat
        condition_names{run} = [SPM.Sess(run).U.name]; %#ok<*AGROW>
        
        % Generate filters for the conditions
        filters(run).(curr_cond_field) = filter_condition(condition_names{run}, curr_cond);
        
        % Get onsets by condition
        curr_ons = sort(vertcat(SPM.Sess(run).U(filters(run).(curr_cond_field)).ons));
        ons(run).(curr_cond_field) = curr_ons;   
        
        % Get causes by condition
        curr_u = sort(vertcat(SPM.Sess(run).U(filters(run).(curr_cond_field)).u));
        u(run).(curr_cond_field) = curr_u;
    end
end

conds = conds'; % transpose

end 