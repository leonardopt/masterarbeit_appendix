function f = our_model(dt, trial_max_val, rest_min_val, randn_noise_factor)

% randn_noise_factor: set it to 0 if you don't want to add any noise to the
% model 

t  = [0:dt:90]; % s
f  = zeros(size(t)); % function to generate
% trial_max_val = 1;
% rest_min_val = -1;

% Initialise figure 
fig = figure('Name', 'Our model');
plot(t, f)

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

plot(t, f)


%% Function 3
% rest1-rest2: 54-56 seconds
t_rest2 = t_rest1 + 2;
filter_3 = (t > t_rest1) & (t <= t_rest2);
curr_t = t(filter_3);
curr_t = curr_t - curr_t(1); % set start to 0 

f(filter_3) = linspace(trial_max_val, rest_min_val, length(curr_t)); % linear from 0 to 1 after cueon s
plot(t, f);


%% Function 4
% rest2-end: 56-90 seconds 
filter_4 = t >= t_rest2;
curr_t = t(filter_4);
curr_t = curr_t - curr_t(1); % set start to 0 

f4 = logspace(0, rest_min_val, length(curr_t));
norm_f4 = (f4 - min(f4) )/ (max(f4) - min(f4));
f(filter_4) = norm_f4 *rest_min_val; % linear from 0 to 1 after cueon s


%% Add random noise if selected        
f = f + randn_noise_factor* randn(size(f));


%% Plot 
plot(t, f);

%% Arrange axis
axis([t(1) t(end) rest_min_val*1.5 trial_max_val*1.5]);
xlabel('dt');
ylabel('modelled activity');
%legend(sprintf('dt= %0.2f', dt), 'Location', 'eastoutside')
title(sprintf('dt= %g', dt))



end 



