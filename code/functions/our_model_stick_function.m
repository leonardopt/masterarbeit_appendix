
% Function to create a model stick function for one block
% codition: cue, trial, rest
% trialmodality: fadeoff or suppression

function f = our_model_stick_function(dt, condition, response, trialendmodality, trailendweight)

if ~exist('trialendmodality', 'var')
    trialmodality = 'none';
end

t  = [0:dt:90]; % s
f  = zeros(size(t)); % function to generate



% Initialise figurea
fig = figure('Name', 'Stick function');

switch condition
    
    case 'cue'
        % Model only the cue
        f(1) = response;
        
        
        
    case 'trial'
        
        % There are 8 trials, one about every 6 seconds
        t_trialon = 2; % s
        t_rest1 = 54; % s
        filter_1 = (t >= t_trialon) & (t <= t_rest1);
        
        filter_ons = find(filter_1 == 1);
        filter_onstr = filter_ons(1);
        filter_onsrs = filter_ons(end);
        
        gap = floor((filter_onsrs - filter_onstr)/ 8); % assumes that trials are presented
        
        ons = [filter_onstr:gap:filter_onsrs];
        ons(end) = '';
        if length(ons) ~= 8
            error('Number of modelled trials is not 8. Please check stick function.')
        end
        
        f(ons) = response;
        
        t_rest1 = 54; % s
        filter_2 = t > t_rest1;
        filter_ons = find(filter_2 == 1);
        filter_onsrs = filter_ons(1);
        filter_endrs = filter_ons(end);
        gap = floor((filter_endrs - filter_onsrs)/ 6);  % according to the causes in the SPM matrices (SPM.Sess.U.u)
        ons = [filter_onsrs:gap:filter_endrs];
        ons(end) = '';
        if length(ons) ~= 6
            error('Number of modelled trials is not 8. Please check stick function.')
        end
        r = repmat(response,size(ons));
        r(end) = 0;     
        
        
        switch trialendmodality
            case 'none'
                r(:) = 0;
            
            case 'fadeoff'
                
                r = r./[1:6]*trailendweight;
                
            case 'suppression'
                
                r = r./[1:6]*(-trailendweight);
                
            otherwise
                error('Unknown trialmodality.')
        end
        f(ons) = r;
        
        
    case 'rest'
        
        t_rest1 = 54; % s
        filter_2 = t > t_rest1;
        filter_ons = find(filter_2 == 1);
        filter_onsrs = filter_ons(1);
        filter_endrs = filter_ons(end);
        
        gap = floor((filter_endrs - filter_onsrs)/ 6);  % according to the causes in the SPM matrices (SPM.Sess.U.u)
        
        ons = [filter_onsrs:gap:filter_endrs];
        ons(end) = '';
        if length(ons) ~= 6
            error('Number of modelled trials is not 8. Please check stick function.')
        end
        
        f(ons) = response;
        
    case 'rest_onestick'
        
        t_rest1 = 54; % s
        filter_2 = t > t_rest1;
        filter_ons = find(filter_2 == 1);
        filter_onsrs = filter_ons(1);
        f(filter_onsrs) = response;
        
    otherwise
        error('Case unknown')
        
        
end

plot(f)

end

