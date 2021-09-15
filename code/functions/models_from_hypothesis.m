% Function to make models from a specified hypothesis
% Hypothesis can be:
% 'boxcar'
% 'sticknoaftereffect_onestickrest'
% 'sticknoaftereffect_allstickrest'
% 'sticknoaftereffect_boxrest'
% 'sticknoaftereffect_randstickrest'
% 'sticknoaftereffect_randstickrest_higherfirststick'
% 'stickfadeoff'
% 'sticksuppression'
% 'stickshortsuppression'
% 'piecewisesuppression'
% 'piecewisefadeoff'

% input:
% hypothesis: string with name of the model to be run
% dt: custom dt. Optional (default 0.0606)
% response: an array with 3 response parameters (cue, trial, rest model).
%           Optional (default is [1 10 10]
% plotresult: plot models. Optional (default is 1)

function [modelcue, modeltrial, modelrest] = models_from_hypothesis(hypothesis,dt, response, plotresult)

%%
% If there's a 'sim-' before the name of the hypothesis, remove it
if strcmp(hypothesis(1:4), 'sim-'), hypothesis = hypothesis(5:end); end
if ~exist('response', 'var'), response = [1 10 10]; end
if ~exist('dt', 'var'), dt = 0.0606; end
if ~exist('plotresult', 'var'), plotresult = 1; end
% if ~exist('bf', 'var'), bf = ''; end


t  = [0:dt:90]; % s
% Unpack response array
% cueresponse_boxcar   = response(1);
% trialresponse_boxcar = response(2);
% restresponse_boxcar  = response(3);

% % Determine height of sticks in relation to boxcar
% freq = 6; % frequency of the sticks (in seconds)
% cueresponse   = cueresponse_boxcar/dt * freq;
% trialresponse = trialresponse_boxcar/dt * freq;
% restresponse  = restresponse_boxcar/dt * freq;

% Unpack response array
cueresponse   = response(1);
trialresponse = response(2);
restresponse  = response(3);

% Determine height of boxcar in relation to sticks
freq = 6; % frequency of the sticks (in seconds)
cueresponse_boxcar   = cueresponse * dt / freq;
trialresponse_boxcar = trialresponse * dt / freq;
restresponse_boxcar  = restresponse * dt / freq;



switch hypothesis
    
    case 'boxcar'
        
        %% Model cue
        %  function to generate
        f  = zeros(size(t));
        % Function 1
        % cue-trial: 0-2 seconds
        t_trialon = 2;
        filter_1 = t < t_trialon;
        f(filter_1) = cueresponse_boxcar;
        % Function 2
        t_trialon = 2;
        filter_2 = (t >= t_trialon);
        f(filter_2) = 0;
        % Store model
        modelcue = f ;
        
        %% Model trial
        %  function to generate
        g  = zeros(size(t));
        % Function 1
        % cue-trial: 0-2 seconds
        t_trialon = 2; % s
        t_rest1 = 54; % s
        gfilter_1 = (t >= t_trialon) & (t <= t_rest1);
        curr_t = t(gfilter_1);
        curr_t = curr_t - curr_t(1); % set start to 0
        g(gfilter_1) = trialresponse_boxcar; % linear from 0 to 1 after cueon s
        % Function 2
        gfilter_2 = t > t_rest1;
        g(gfilter_2) = 0;
        % Store model
        modeltrial = g ;
        
        %% Model rest
        %  function to generate
        h  = zeros(size(t));
        % Function 1
        % cue-trial: 0-2 seconds
        t_cueon = 0; % s
        t_rest1 = 54; % s
        hfilter_1 = (t > t_cueon) & (t <= t_rest1);
        curr_t = t(hfilter_1);
        curr_t = curr_t - curr_t(1); % set start to 0
        h(hfilter_1) = 0; % linear from 0 to 1 after cueon s
        % Function 2
        hfilter_2 = t > t_rest1;
        h(hfilter_2) = restresponse_boxcar;
        % Store mode
        modelrest = h;
        
        
    case 'sticknoaftereffect_onestickrest'
        
        f  = zeros(size(t)); % function to generate
        
        %% Model cue
        modelcue = f;
        modelcue(1) = cueresponse;
        %% Model trial
        % There are 8 trials, one about every 6 seconds
        f  = zeros(size(t));
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
        
        f(ons) = trialresponse;
        
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
        r = repmat(trialresponse,size(ons));
        r(end) = 0;
        % Trial end
        r(:) = 0;
        f(ons) = r;
        modeltrial = f;
        
        %% Model rest
        f  = zeros(size(t));
        t_rest1 = 54; % s
        filter_2 = t > t_rest1;
        filter_ons = find(filter_2 == 1);
        filter_onsrs = filter_ons(1);
        f(filter_onsrs) = restresponse;
        modelrest = f;
        
    case 'sticknoaftereffect_allstickrest'
        f  = zeros(size(t)); % function to generate
        
        %% Model cue
        modelcue = f;
        modelcue(1) = cueresponse;
        %% Model trial
        % There are 8 trials, one about every 6 seconds
        f  = zeros(size(t));
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
        
        f(ons) = trialresponse;
        
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
        r = repmat(trialresponse,size(ons));
        r(end) = 0;
        % Trial end
        r(:) = 0;
        f(ons) = r;
        modeltrial = f;
        
        %% Model rest
        f  = zeros(size(t));
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
        
        f(ons) = restresponse;
        modelrest = f;
        
    case 'sticknoaftereffect_randstickrest'
        
        f  = zeros(size(t)); % function to generate
        
        %% Model cue
        modelcue = f;
        modelcue(1) = cueresponse;
        %% Model trial
        % There are 8 trials, one about every 6 seconds
        f  = zeros(size(t));
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
        
        f(ons) = trialresponse;
        
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
        r = repmat(trialresponse,size(ons));
        r(end) = 0;
        % Trial end
        r(:) = 0;
        f(ons) = r;
        modeltrial = f;
        
        %% Model rest
        
        f  = zeros(size(t));
        t_rest1 = 54; % s
        filter_2 = t > t_rest1;
        filter_ons = find(filter_2 == 1);
        filter_onsrs = filter_ons(1);
        filter_endrs = filter_ons(end);
        
        % Number ofrandom sticks
        n_sticks = 4;
        % Stick probabiliy
        x = 1:size(t(filter_onsrs:end), 2);
        y = 1./x;
        % Cum dist function
        d = cumsum(y)/sum(y);
        
        % Sticks should be at least n seconds from each other
        % get how many elements of t correspond to n seconds
        minsec = 6; % frequency of trials
        mindistance = numel(0:dt:minsec);
        
        while 1
            % Sample four sticks randomly from distribution
            sticks = rand(n_sticks, 1);
            [minval, closestidx] = min(abs(d - sticks), [], 2);
            closestval = d(closestidx); % to avoid mistakes due to floats
            % Find time points given the sampled probabilities
            time_sticks = find(ismember(d, closestval));
            % Conditions: numbers sampled were actually found in the array
            % AND the distance between the time points is bigger than the
            % minimum distance of n seconds
            if length(time_sticks)== n_sticks && sum(diff(time_sticks)<mindistance)==0
                break
            end
        end
        
        % These are the onsets of our rest model sticks
        ons = filter_onsrs + time_sticks;
        
        % Store results in model
        f(ons) = restresponse;
        modelrest = f;
        
    case 'sticknoaftereffect_randstickrest_higherfirststick'
        
        f  = zeros(size(t)); % function to generate
        
        %% Model cue
        modelcue = f;
        modelcue(1) = cueresponse;
        %% Model trial
        % There are 8 trials, one about every 6 seconds
        f  = zeros(size(t));
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
        
        f(ons) = trialresponse;
        
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
        r = repmat(trialresponse,size(ons));
        r(end) = 0;
        % Trial end
        r(:) = 0;
        f(ons) = r;
        
        % Make first stick 50% higher 
        sticks = find(f>0);
        f(sticks(1)) = trialresponse + trialresponse/2;
        
        
        modeltrial = f;
        
        %% Model rest
        
        f  = zeros(size(t));
        t_rest1 = 54; % s
        filter_2 = t > t_rest1;
        filter_ons = find(filter_2 == 1);
        filter_onsrs = filter_ons(1);
        filter_endrs = filter_ons(end);
        
        % Number ofrandom sticks
        n_sticks = 4;
        % Stick probabiliy
        x = 1:size(t(filter_onsrs:end), 2);
        y = 1./x;
        % Cum dist function
        d = cumsum(y)/sum(y);
        
        % Sticks should be at least n seconds from each other
        % get how many elements of t correspond to n seconds
        minsec = 6; % frequency of trials
        mindistance = numel(0:dt:minsec);
        
        while 1
            % Sample four sticks randomly from distribution
            sticks = rand(n_sticks, 1);
            [minval, closestidx] = min(abs(d - sticks), [], 2);
            closestval = d(closestidx); % to avoid mistakes due to floats
            % Find time points given the sampled probabilities
            time_sticks = find(ismember(d, closestval));
            % Conditions: numbers sampled were actually found in the array
            % AND the distance between the time points is bigger than the
            % minimum distance of n seconds
            if length(time_sticks)== n_sticks && sum(diff(time_sticks)<mindistance)==0
                break
            end
        end
        
        % These are the onsets of our rest model sticks
        ons = filter_onsrs + time_sticks;
        
        % Store results in model
        f(ons) = restresponse;
        modelrest = f;
        
        
    case 'sticknoaftereffect_boxrest'
        
        f  = zeros(size(t)); % function to generate
        
        %% Model cue
        modelcue = f;
        modelcue(1) = cueresponse;
        
        %% Model trial
        % There are 8 trials, one about every 6 seconds
        f  = zeros(size(t));
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
        
        f(ons) = trialresponse;
        
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
        r = repmat(trialresponse,size(ons));
        r(end) = 0;
        % Trial end
        r(:) = 0;
        f(ons) = r;
        modeltrial = f;
        
        %% Model rest
        %  function to generate
        h  = zeros(size(t));
        % Function 1
        % cue-trial: 0-2 seconds
        t_cueon = 0; % s
        t_rest1 = 54; % s
        hfilter_1 = (t > t_cueon) & (t <= t_rest1);
        curr_t = t(hfilter_1);
        curr_t = curr_t - curr_t(1); % set start to 0
        h(hfilter_1) = 0; % linear from 0 to 1 after cueon s
        % Function 2
        hfilter_2 = t > t_rest1;
        h(hfilter_2) = restresponse_boxcar;
        % Store mode
        modelrest = h;
        
        
        
    case 'sticknoaftereffect_boxtrialrest'
        
        f  = zeros(size(t)); % function to generate
        
        %% Model cue
        modelcue = f;
        modelcue(1) = cueresponse;
        
        %% Model trial
        % There are 8 trials, one about every 6 seconds
        f  = zeros(size(t));
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
        
        f(ons) = trialresponse;
        
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
        r = repmat(trialresponse,size(ons));
        r(end) = 0;
        % Trial end
        r(:) = 0;
        f(ons) = r;
        modeltrial = f;
        
        %% Model rest
        %  function to generate
        h  = zeros(size(t));
        h(:) = restresponse_boxcar;
        modelrest = h;
        
        case 'sticknoaftereffect_boxallblock'
        
        f  = zeros(size(t)); % function to generate
        
        %% Model all block 
        
        f(:)= restresponse_boxcar;    
        modelcue = f;
        
        %% Model trial
        % There are 8 trials, one about every 6 seconds
        f  = zeros(size(t));
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
        
        f(ons) = trialresponse;
        
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
        r = repmat(trialresponse,size(ons));
        r(end) = 0;
        % Trial end
        r(:) = 0;
        f(ons) = r;
        modeltrial = f;
        
        %% Model rest
        %  function to generate
        h  = zeros(size(t));
        % Function 1
        % cue-trial: 0-2 seconds
        t_cueon = 0; % s
        t_rest1 = 54; % s
        hfilter_1 = (t > t_cueon) & (t <= t_rest1);
        curr_t = t(hfilter_1);
        curr_t = curr_t - curr_t(1); % set start to 0
        h(hfilter_1) = 0; % linear from 0 to 1 after cueon s
        % Function 2
        hfilter_2 = t > t_rest1;
        h(hfilter_2) = restresponse_boxcar;
        % Store mode
        modelrest = h;
        
    case 'stickfadeoff'
        
        
        f  = zeros(size(t)); % function to generate
        
        %% Model cue
        modelcue = f;
        modelcue(1) = cueresponse;
        %% Model trial
        % There are 8 trials, one about every 6 seconds
        f  = zeros(size(t));
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
        
        f(ons) = trialresponse;
        
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
        r = repmat(trialresponse,size(ons));
        r(end) = 0;
        % Trial end
        r = r./[1:6];
        f(ons) = r;
        modeltrial = f;
        
        %% Model rest
        
        f  = zeros(size(t));
        t_rest1 = 54; % s
        filter_2 = t > t_rest1;
        filter_ons = find(filter_2 == 1);
        filter_onsrs = filter_ons(1);
        filter_endrs = filter_ons(end);
        
        % Number ofrandom sticks
        n_sticks = 4;
        % Stick probabiliy
        x = 1:size(t(filter_onsrs:end), 2);
        y = 1./x;
        % Cum dist function
        d = cumsum(y)/sum(y);
        
        % Sticks should be at least n seconds from each other
        % get how many elements of t correspond to n seconds
        minsec = 6; % frequency of trials
        mindistance = numel(0:dt:minsec);
        
        while 1
            % Sample four sticks randomly from distribution
            sticks = rand(n_sticks, 1);
            [minval, closestidx] = min(abs(d - sticks), [], 2);
            closestval = d(closestidx); % to avoid mistakes due to floats
            % Find time points given the sampled probabilities
            time_sticks = find(ismember(d, closestval));
            % Conditions: numbers sampled were actually found in the array
            % AND the distance between the time points is bigger than the
            % minimum distance of n seconds
            if length(time_sticks)== n_sticks && sum(diff(time_sticks)<mindistance)==0
                break
            end
        end
        
        % These are the onsets of our rest model sticks
        ons = filter_onsrs + time_sticks;
        
        % Store results in model
        f(ons) = restresponse;
        modelrest = f;
        
        
        
    case 'sticksuppression'
        
        f  = zeros(size(t)); % function to generate
        
        %% Model cue
        modelcue = f;
        modelcue(1) = cueresponse;
        %% Model trial
        % There are 8 trials, one about every 6 seconds
        f  = zeros(size(t));
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
        
        f(ons) = trialresponse;
        
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
        r = repmat(trialresponse,size(ons));
        r(end) = 0;
        % Trial end
        r = -r./[1:6];
        
        f(ons) = r;
        modeltrial = f;
        
        %% Model rest
        %% Model rest
        
        f  = zeros(size(t));
        t_rest1 = 54; % s
        filter_2 = t > t_rest1;
        filter_ons = find(filter_2 == 1);
        filter_onsrs = filter_ons(1);
        filter_endrs = filter_ons(end);
        
        % Number ofrandom sticks
        n_sticks = 4;
        % Stick probabiliy
        x = 1:size(t(filter_onsrs:end), 2);
        y = 1./x;
        % Cum dist function
        d = cumsum(y)/sum(y);
        
        % Sticks should be at least n seconds from each other
        % get how many elements of t correspond to n seconds
        minsec = 6; % frequency of trials
        mindistance = numel(0:dt:minsec);
        
        while 1
            % Sample four sticks randomly from distribution
            sticks = rand(n_sticks, 1);
            [minval, closestidx] = min(abs(d - sticks), [], 2);
            closestval = d(closestidx); % to avoid mistakes due to floats
            % Find time points given the sampled probabilities
            time_sticks = find(ismember(d, closestval));
            % Conditions: numbers sampled were actually found in the array
            % AND the distance between the time points is bigger than the
            % minimum distance of n seconds
            if length(time_sticks)== n_sticks && sum(diff(time_sticks)<mindistance)==0
                break
            end
        end
        
        % These are the onsets of our rest model sticks
        ons = filter_onsrs + time_sticks;
        
        % Store results in model
        f(ons) = restresponse;
        modelrest = f;
        
        
    case 'stickshortsuppression'
        
        f  = zeros(size(t)); % function to generate
        
        %% Model cue
        modelcue = f;
        modelcue(1) = cueresponse;
        %% Model trial
        % There are 8 trials, one about every 6 seconds
        f  = zeros(size(t));
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
        
        f(ons) = trialresponse;
        
        t_rest1 = 54; % s
        filter_2 = t > t_rest1;
        filter_ons = find(filter_2 == 1);
        filter_ons = filter_ons(1);
        
        f(filter_ons) = -trialresponse;
        
        modeltrial = f;
        
        %% Model rest
        %% Model rest
        
        f  = zeros(size(t));
        t_rest1 = 54; % s
        filter_2 = t > t_rest1;
        filter_ons = find(filter_2 == 1);
        filter_onsrs = filter_ons(1);
        filter_endrs = filter_ons(end);
        
        % Number ofrandom sticks
        n_sticks = 4;
        % Stick probabiliy
        x = 1:size(t(filter_onsrs:end), 2);
        y = 1./x;
        % Cum dist function
        d = cumsum(y)/sum(y);
        
        % Sticks should be at least n seconds from each other
        % get how many elements of t correspond to n seconds
        minsec = 6; % frequency of trials
        mindistance = numel(0:dt:minsec);
        
        while 1
            % Sample four sticks randomly from distribution
            sticks = rand(n_sticks, 1);
            [minval, closestidx] = min(abs(d - sticks), [], 2);
            closestval = d(closestidx); % to avoid mistakes due to floats
            % Find time points given the sampled probabilities
            time_sticks = find(ismember(d, closestval));
            % Conditions: numbers sampled were actually found in the array
            % AND the distance between the time points is bigger than the
            % minimum distance of n seconds
            if length(time_sticks)== n_sticks && sum(diff(time_sticks)<mindistance)==0
                break
            end
        end
        
        % These are the onsets of our rest model sticks
        ons = filter_onsrs + time_sticks;
        
        % Store results in model
        f(ons) = restresponse;
        modelrest = f;
        
    case 'piecewisesuppression'
        
        f  = zeros(size(t)); % function to generate
        trial_max_val = trialresponse_boxcar;
        rest_min_val  = -trialresponse_boxcar;
        
        
        %% Function 1
        % cue-trial: 0-2 seconds
        t_cueon = 0; % s
        t_trialon = 2; % s
        filter_1 = (t > t_cueon) & (t <= t_trialon);
        curr_t = t(filter_1);
        curr_t = curr_t - curr_t(1); % set start to 0
        f(filter_1) = linspace(0, trial_max_val, length(curr_t)); % linear from 0 to 1 after cueon s
        
        plot(t, f)
        
        %% Function 2
        % trial-rest: 2-54 seconds
        t_rest1 = 54; % s
        filter_2 = (t > t_trialon) & (t <= t_rest1);
        f(filter_2) = trial_max_val + 0 * t(filter_2); % linear from 0 to 1 after cueon s
        
        
        
        %% Function 3
        % rest1-rest2: 54-56 seconds
        t_rest2 = t_rest1 + 2;
        filter_3 = (t > t_rest1) & (t <= t_rest2);
        curr_t = t(filter_3);
        curr_t = curr_t - curr_t(1); % set start to 0
        
        f(filter_3) = linspace(trial_max_val, rest_min_val, length(curr_t)); % linear from 0 to 1 after cueon s
        
        
        
        %% Function 4
        % rest2-end: 56-90 seconds
        filter_4 = t >= t_rest2;
        curr_t = t(filter_4);
        curr_t = curr_t - curr_t(1); % set start to 0
        
        f4 = logspace(0, rest_min_val, length(curr_t));
        norm_f4 = (f4 - min(f4) )/ (max(f4) - min(f4));
        f(filter_4) = norm_f4 *rest_min_val; % linear from 0 to 1 after cueon s
        
        modelcue   = f;
        modeltrial = f;
        modelrest  = f;
        
    case 'piecewisefadeoff'
        
        % function to generate
        f  = zeros(size(t));
        trial_max_val = trialresponse_boxcar;
        rest_min_val = 0;
        
        %% Function 1
        % cue-trial: 0-2 seconds
        t_cueon = 0; % s
        t_trialon = 2; % s
        filter_1 = (t > t_cueon) & (t <= t_trialon);
        curr_t = t(filter_1);
        curr_t = curr_t - curr_t(1);
        f(filter_1) = linspace(0, trial_max_val, length(curr_t)); % linear from 0 to 1 after cueon s
        
        %% Function 2
        % trial-rest: 2-54 seconds
        t_rest1 = 54; % s
        filter_2 = (t > t_trialon) & (t <= t_rest1);
        f(filter_2) = trial_max_val + 0 * t(filter_2);
        
        
        
        %% Function 3
        % rest1-rest2: 54-56 seconds
        t_rest2 = t_rest1 + t(end);
        filter_3 = (t > t_rest1) & (t <= t_rest2);
        curr_t = t(filter_3);
        curr_t = curr_t - curr_t(1);
        
        f(filter_3) = linspace(trial_max_val, rest_min_val, length(curr_t));
        
        modelcue   = f;
        modeltrial = f;
        modelrest  = f;
        
    otherwise
        allmodelnames = {'boxcar', ...
            'sticknoaftereffect_onestickrest', ...
            'sticknoaftereffect_allstickrest', ...
            'sticknoaftereffect_boxrest', ...
            'stickfadeoff', ...
            'sticksuppression', ...
            'stickshortsuppression', ...
            'piecewisesuppression', ...
            'piecewisefadeoff'};
        
        warning('Invalid model name')
        fprintf('Valid model names: % s\n', allmodelnames{:})      
        error('Model name unknown')
        
        
end


%
% %% Plot results
%
% if plotresult
%     subplot(2,2,1)
%     plot(modelcue)
%
%     subplot(2,2,2)
%     plot(modeltrial)
%
%     subplot(2,2,3)
%     plot(modelrest)
%
%
% end









end