% Gets onsets, causes, etc. but grouped by condition
%
% Kai/Leonardo ~2019

function [U, filters] = group_U_conditions(SPM, run)

% All conditions 
conds = {'left', 'right', 'aeh-s', 'aeh-c', 'kon-s', 'kon-c', 'sum-s', 'sum-c'};
filters = struct();

% Get condition names by run
for i = 1:length(conds)
    
    % Get current condition
    curr_cond = conds{i};
    curr_cond_field = strrep(curr_cond, '-', '_');
    
    % Extract names of the condition from the SPM mat
    condition_names{run} = [SPM.Sess(run).U.name];
    
    % Generate filters for the conditions
    filter = filter_condition(condition_names{run}, curr_cond);
    
    % Get onsets by condition
    curr_ons = sort(vertcat(SPM.Sess(run).U(filter).ons));
    
    % Get causes by condition
    curr_u = sort(vertcat(SPM.Sess(run).U(filter).u));
    
    u.(curr_cond_field) = curr_u;
     
    U(i).name = {curr_cond};
    U(i).ons  = curr_ons;
    U(i).u    = curr_u;
    U(i).dur  = zeros(size(curr_ons));
    U(i).orth = 1;
    U(i).P    = struct('name', 'none', 'h', 0, 'i', 1);
    U(i).dt   = SPM.Sess(1).U(1).dt;
    
    % Store filters 
    filters(run).(curr_cond_field) = filter;
end