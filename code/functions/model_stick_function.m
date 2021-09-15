% Function to create a model stick function for one block.
%
% IN
%   cfg with:
%    folder path to save figure to (e.g. run_cfg with "sim-") folder name
%    cfg.saa.data_gen.trialboxcar = 1/0 (default = 0, set to sticks)
%    cfg.saa.data_gen.restboxcar = 1/0 OR
%    cfg.saa.data_gen.resttype = 'reststicks', 'restboxcar', 'randsticks', 
%     'lowfreq'
%    cfg.saa.data_gen.cuttrialend    = 1;
%    cfg.saa.data_gen.trialendstick  = 1; % applies to some aftereffects 
%      to mimic HRF effect (prevents sudden drop after last trial onset 
%      within task phase)
%    cfg.saa.data_gen.trialoffset    = 2; % in sec. not TR
%    cfg.saa.data_gen.modelrest      = 1; % <--- CHECK
%    cfg.saa.data_gen.negativerest   = 0; % positive/negative rest signal
%
%   trialendmodality ('fadeoff', 'suppression', 'shortsuppression', 'none')
%    no default
%
% Optional IN
%   dt: TR/n_slices. Default (for afx) = 0.0606;
%      task+rest duration = 90 seconds, TR = 2; 45 timepoints per block;
%      276 images total per run
%   trialendweight: default = 1; linear scaling of endpoints (0 - Inf)
%   response: signal amplitude  parameter? vector of three values for
%      cue, trial and rest in that order;
%      e.g. [1 10 10] = 1 for cue, 10 for trial and rest
%
% OUT
%   model parameters of simulated timecourse for cue, task and rest
%
% Created by Leonardo, ~2019
% Adapted by Ingmar, 06/2020

function [modelcue, modeltrial, modelrest] = model_stick_function(cfg, trialendmodality, trialendweight, response, dt)
%% Check input

% No weighting for after trial (e.g. for larger fadeoff/suppression)
if ~exist('trialendweight', 'var'), trialendweight = 1; end

% Standard response is equal for task and rest
if ~exist('response', 'var'), response = [1 10 10]; end

% TR / slices (for slice timing correction)
if ~exist('dt', 'var'), dt = 2/33; end

% Checks
assert(logical(exist('cfg', 'var')), ...
    'Saving stick function plot will fail without cfg with simulation directory.');
assert(isfield(cfg, 'saa'), ...
    'cfg.saa field (and cfg.saa.data_gen) are required.');
assert(logical(exist('trialendmodality', 'var')), ...
    'Please enter trialendmodality as first argument (e.g. "fadeoff")');

%% Other settings to default to

% -------------------% Task trials
if isfield(cfg.saa.data_gen, 'trialboxcar') && cfg.saa.data_gen.trialboxcar
    trialendmodality = 'none'; % Not meant for aftereffect simulation!
else
    cfg.saa.data_gen.trialboxcar = 0;
end

% -------------------% Rest
if ~isfield(cfg.saa.data_gen, 'modelrest')
    cfg.saa.data_gen.modelrest = 1; 
end

if ~isfield(cfg.saa.data_gen, 'negativerest')
    cfg.saa.data_gen.negativerest = 0;
end

if isfield(cfg.saa.data_gen, 'restboxcar') && cfg.saa.data_gen.restboxcar
    cfg.saa.data_gen.resttype = 'boxcar'; % compatibility
end

% -------------------% Plotting
if ~isfield(cfg.saa.data_gen, 'plotting')
    cfg.saa.data_gen.plotting = 0;
end

disp('    Simulation model settings:'); disp(cfg.saa.data_gen);

%% Onsets & filters

% Specific to AFX, do not change ------- %
time                = 0:dt:90; % timecourse
func                = zeros(size(time)); % preallocate
t_trialonset        = 2; % 1 TR
t_trialend          = 54; % 27 TR
t_restonset         = 54; % sec
t_trial_unadj       = t_trialonset;
t_trialend_unadj    = t_trialend;
% -------------------------------------- %

if isfield(cfg.saa.data_gen, 'trialoffset') && ~isempty(cfg.saa.data_gen.trialoffset)
    trialoffset     = cfg.saa.data_gen.trialoffset; % in sec
    t_trialonset    = t_trialonset + trialoffset ;
    t_trialend      = t_trialend + trialoffset;
    fprintf('    Modeltrial is offset with %i second(s)! Might better mimic empirical results.\n', trialoffset);
end

% Filters
trial_filter        = (time >= t_trialonset) & (time <= t_trialend);
trial_tp_filter     = find(trial_filter == 1);
rest_filter         = time > t_restonset;
rest_tp_filter      = find(rest_filter == 1);

% Keep aftereffect timing the same no matter any offset, except if adding
% extra stick
t_filter_adj        = (time >= t_trial_unadj) & (time <= t_trialend_unadj);
t_tp_filter_adj     = find(t_filter_adj == 1);

%% Cue

modelcue    = func;
modelcue(1) = response(1);

%% Trial phase
% There are 8 trials, start every 6.5 seconds (ITI)

% Sticks at trial onsets
trial_stick_gap         = floor((rest_tp_filter(1) - trial_tp_filter(1))/8);
trial_stick             = trial_tp_filter(1):trial_stick_gap:rest_tp_filter(1);
modeltrial              = func;
trial_stick(end)        = '';
modeltrial(trial_stick) = response(2);

% Check
assert(length(trial_stick) == 8, ...
    'Number of modelled trials is not 8. Check stick function.');

% Boxcar
if cfg.saa.data_gen.trialboxcar
    modeltrial(trial_filter) = response(2)/trial_stick_gap; % boxcar
    fprintf('    Trial boxcar level was reduced to 1/%i of sticks\n', trial_stick_gap);
end

%% Trial end phase

% New code to keep rest "pulses" equal to task onsets
trialend_stick = t_tp_filter_adj(end):trial_stick_gap:rest_tp_filter(end);

% Trial end modality
fprintf('    Simulation: trialendweight of trialendmodality "%s" = %02f\n\n', ...
    trialendmodality, trialendweight)

switch trialendmodality
    
    case 'none'
        
        fprintf('    No aftereffect selected ("none")...\n');
        
    case 'linearfadeoff'
        
        x = (response(2)/length(trialend_stick));
        y = (response(2)-x):-x:0;
        modeltrial(trialend_stick) = y;
        
    case 'fadeoff'
        % defined now as trial level log-asymptotic (weighted) to 1
        
        logstick = 1; % for calculating base log, 1 stick extra
        
        log_fadeoff = logspace(log10(response(2)), log10(1), ...
            length(trialend_stick) + logstick) - 1;
        
        % cut first response, which is same level as trial
        fadeoff = log_fadeoff(2:end);
        
        % add weight factor
        fadeoff = fadeoff .* trialendweight;
        
        if cfg.saa.data_gen.cuttrialend && length(fadeoff) > 5 % cut sticks to 5
            fprintf('    Limiting aftereffect sticks to 5\n');
            fadeoff         = fadeoff(1:5);
            trialend_stick  = trialend_stick(1:5);
        end
        
        modeltrial(trialend_stick) = fadeoff;
        
    case {'sustainedsuppression', 'suppression', 'inhibition'}
        % negative trial level log-asymptotic (weighted) to zero
        
        logstick = 0; % for calculating base log, 1 stick extra
        
        % log10 below 1 is "nonlinear", so add 1 after inversion
        log_suppression = -logspace(log10(response(2) + 1), ...
            log10(1) .* trialendweight, length(trialend_stick) + logstick) + 1;
        
        if cfg.saa.data_gen.cuttrialend && length(log_suppression) >= length(trialend_stick)
            fprintf('    Limiting aftereffect sticks to %i\n', length(trialend_stick));
            if cfg.saa.data_gen.trialendstick % mimics HRF fadeout to baseline before sign flip of suppr.
                suppression     = log_suppression;
                suppression(1)  = 0;
            else
                suppression = log_suppression(2:end); % cut first response
            end
            suppression     = suppression(1:5);
            trialend_stick  = trialend_stick(1:5);
        end
        
        modeltrial(trialend_stick) = suppression;
        
    case 'shortsuppression'
        % Quick return to baseline asymptote (log10 below 1 is
        % "nonlinear", so add 1 after inversion)
        
        logstick = 1; % for calculating base log, 1 stick extra
        
        x = 2; % 2 sticks
        log_suppression = -logspace(log10(response(2) + 1), ...
            log10(1) .* trialendweight, x + logstick) + 1; % length + 1 to end at 0
        shortsuppression     = log_suppression(1:2); % select 2 sticks
        shortsuppression     = shortsuppression .* trialendweight; % add weighting factor
        
        if cfg.saa.data_gen.trialendstick
            modeltrial(trialend_stick(1))   = response(2)/2;
            modeltrial(trialend_stick(2:3)) = shortsuppression;
        else
            modeltrial(trialend_stick(1:2)) = shortsuppression;
        end
        
    case 'gradualsuppression' % Gradual drop below baseline and earlier return with 3rd order poly (z-curve)
        
        x = 1:4; % sticks
        gradsuppr = (1.3.*x.^3) - (5.8 .* x.^2) + 10;
        gradsuppr = gradsuppr./(10/response(2)); % calibrated on y = 10
        modeltrial(trialend_stick(1:4)) = gradsuppr;
        
        % 2nd order polynomial (u-curve)
        % sticklen    = 0:4;
        % gradsuppr   = (3 .* sticklen.^2) - (12 .* sticklen) + response(2);
        
    otherwise
        error('%s is an unknown trial end modality.', trialendmodality)
end

disp('    Unique trial model stick values:');
disp(fliplr(unique(modeltrial)));

%% Rest model (same timepoints as post-trial)

% Default = 1
if ~isfield(cfg.saa.data_gen, 'modelrest') || cfg.saa.data_gen.modelrest 
    restresponse = response(3);
        
    if restresponse == 0 % Still no model
        modelrest(rest_filter) = 0;
        disp('    cfg.saa.data_gen.modelrest == 1 but model response == 0');
    end
    
    switch cfg.saa.data_gen.resttype
        case 'reststicks' 
            % Sticks sampled at same rate as trials
            modelrest(trialend_stick) = restresponse;
            
        case 'restboxcar' 
            % Boxcar model (single sustained response) spanning total rest length
            modelrest(rest_filter) = response(3)/trial_stick_gap; % boxcar
            fprintf('    Stick function: rest boxcar level reduced to 1/%ith of input\n', trial_stick_gap);
            
        case 'randsticks' 
            % A few sticks randomized, custom distribution (hyperbolic?)
            
            nsticks     = length(trial_stick);
            restlength  = length(rest_tp_filter);
            distrib     = rest_tp_filter + (restlength-rest_tp_filter)/restlength; % distrib(1:16)  = abs(8-1:16) % other possible distrib
            dist_cdf    = cumsum(distrib)/sum(distrib); % cumulative sum
            uni_sample  = rand(nsticks, 1);
            x_sample    = nan(size(uni_sample));
            for uni_ind = 1:length(uni_sample)
                x_sample(uni_ind) = find(~(dist_cdf <= uni_sample(uni_ind)), 1);
            end
            x_sample = x_sample + rest_tp_filter(1);
            modelrest(x_sample) = restresponse;
            disp('    Rest sticks randomly sampled (only makes sense to randomize for every "subject"!). Values:'); disp(x_sample);
            
        case 'lowfreq' 
            % Sticks sampled at eg. 0.08 Hz, related to resting state signals in
            % fMRI, onset at peak (half frequency) after rest start

            hz          = 0.08; % randn(1)* 0.01 + 0.04;
            restlength  = length(rest_tp_filter);
            sticksfreq  = 0:(1/(dt*hz)):restlength; % given freq
            onset       = mean(sticksfreq(1:2));
            sticks      = onset:(1/(dt*hz)):restlength; % given freq
            sticks      = round(sticks + rest_tp_filter(1));
            modelrest(sticks) = restresponse;
            fprintf('    Low-frequency rest phase sticks sampled at %g Hz\n\n', hz);
            
            % sd      = sticksfreq(2)/4;
            % onset   = round((randn(1) * sd + (sticksfreq(2)/2))); % randomize
    end
end

% Below baseline
if isfield(cfg.saa.data_gen, 'negativerest') && cfg.saa.data_gen.negativerest
    modelrest = modelrest .* -1;
end

%% Plot sticks

if cfg.saa.data_gen.plotting
    stickfig = figure('name', 'Stick function model', 'Position', [400 350 400 450]);
    
    % Prepare TS
    hold on;
    plot_modelcue = modelcue;
    plot(plot_modelcue, 'g-');
    
    if cfg.saa.data_gen.modelrest
        plot_modelrest = modelrest;
        plot_modelrest(end + 1) = 0;
        restplotstr = 'b-';
        if strcmpi(cfg.saa.data_gen.resttype, 'boxcar')
            plot_modelrest = plot_modelrest .* trial_stick_gap;
            restplotstr = 'b--';
            fprintf('Stick function plot: rest model visibly scaled up to conv. height.\n')
        end
        plot(plot_modelrest, restplotstr);
    end
    
    plot(modeltrial, 'r-'); % plot
    
    % set ticks and limits
    xticks(round(linspace(0, rest_tp_filter(end), 6)));
    ylims = get(gca, 'Ylim');
    set(gca, 'Ylim', [ylims(1)*1.2 ylims(2)*1.4]); % extend figure bounds
    set(gca, 'Xlim', [-25, Inf]);
    
    % set labels
    ylabel('Modelled activity (A.U.)'); xlabel('dt (TR / slices)');
    title('Model spike intervals (stick function)');
    legend('Cue', 'Rest', 'Trial', 'FontSize', 8, 'Orientation', 'horizontal');
    hold off;
    
    %% Save fig
    
    try
        targetfile = fullfile(cfg.derivative_dir, cfg.sim_name, ...
            ['plot-sim_stickfunction-' trialendmodality]);
        if ~isfolder(fullfile(cfg.derivative_dir, cfg.sim_name))
            mkdir(fullfile(cfg.derivative_dir, cfg.sim_name))
        end
        disp('Saving stick function plot to: '); disp(targetfile);
        
        if isfield(cfg.saa.data_gen, 'save_svg') && cfg.saa.data_gen.save_svg
            save_plot(targetfile, stickfig)
        else
            % no save_fig to save memory
            saveas(stickfig, [targetfile '.png']); % Add .png otherwise saved as .fig
        end
    catch e
        warning('Saving plot failed, see error...')
        disp(e.message);
    end
    
end