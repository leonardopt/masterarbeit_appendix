% Function to plot mean vs CI95


function afx_tool_plot_mean_vs_CI95(cfg_sim, plot_name, all_means, all_CI95, odd_labelnames, even_labelnames )

set(0, 'DefaultFigureVisible', 'off')

targetdir                        = cfg_sim.general_info.target_directory;
save_figures                     = cfg_sim.general_info.save_figures;
data_type                        = cfg_sim.general_info.data_type;
voxelrange                       = cfg_sim.general_info.voxelrange;
mean_n_voxels                    = cfg_sim.general_info.mean_n_voxels;


set_figure(cfg_sim.general_info, plot_name, [868   260   436   779])
% Create mask for colors
odd_labelnames = odd_labelnames' ;
even_labelnames = even_labelnames';

filter_odd_trials = ~cellfun(@isempty , strfind(odd_labelnames, 'Trial'));
filter_even_trials = ~cellfun(@isempty , strfind(even_labelnames, 'Trial'));
filter_odd_rest = ~cellfun(@isempty , strfind(odd_labelnames, 'Rest'));
filter_even_rest = ~cellfun(@isempty , strfind(even_labelnames, 'Rest'));

% Get filters for different phases (trial-trial, rest-rest, trial-rest)
all_means_filt_rest  = zeros(size(all_means));
all_means_filt_trial = zeros(size(all_means));
all_means_filt_trialrest = zeros(size(all_means));
for i=1:size(all_means,1)
    for i2 = 1:size(all_means,2)
        all_means_filt_rest(i, i2)  = and(filter_odd_rest(i), filter_even_rest(i2));
        all_means_filt_trial(i, i2) = and(filter_odd_trials(i), filter_even_trials(i2));
        all_means_filt_trialrest(i, i2) = or(and(filter_odd_trials(i), filter_even_rest(i2)),and(filter_odd_rest(i), filter_even_trials(i2))) ;
        
    end
end

all_means_groups                              =  ones(size(all_means)); % cues will be part of trials 
all_means_groups(all_means_filt_trial==1)     = 1;
all_means_groups(all_means_filt_rest==1)      = 2;
all_means_groups(all_means_filt_trialrest==1) = 3;

%
subplot(2,1,1)
x = all_means(:);
y = all_CI95(:);
groups = all_means_groups(:);


% Plot
plot_mean_CI = gscatter(x, y, groups);

% Plot specifications
colours = [ [102,194,165]; ...
            [141, 160, 203];
            [252, 141,98]...
            ];

colours = colours/255;

% plot_mean_CI(1).DisplayName = 'trial-trial';
% plot_mean_CI(2).DisplayName = 'rest-rest';
% plot_mean_CI(3).DisplayName = 'trial-rest';
% plot_mean_CI(1).MarkerSize = sz;
% plot_mean_CI(2).MarkerSize = sz;
% plot_mean_CI(3).MarkerSize = sz;
% plot_mean_CI(1).Color = colours(1, :);
% plot_mean_CI(2).Color = colours(2, :);
% plot_mean_CI(3).Color = colours(3, :);

plotnames = {'trial-trial', 'rest-rest', 'trial-rest'};
sz = 2;

% Set transparency
for i = 1:length(plot_mean_CI)
    set(plot_mean_CI(i), 'DisplayName', plotnames{i});
    set(plot_mean_CI(i), 'Marker','o', 'MarkerSize',sz, 'MarkerEdgeColor','none', 'MarkerFaceColor',colours(i, :));
    set(plot_mean_CI(i).MarkerHandle, 'FaceColorType','truecoloralpha', 'FaceColorData',uint8( 255* [colours(i, 1);colours(i, 2);colours(i, 3);0.3] ) );
end

%
% Mean trial-trial
m_plot_trial = tanh(mean(atanh(all_means(all_means_groups==1))));
CI_plot_trial = tanh(mean(atanh(all_CI95(all_means_groups==1))));
% Mean rest-rest
m_plot_rest = tanh(mean(atanh(all_means(all_means_groups==2))));
CI_plot_rest = tanh(mean(atanh(all_CI95(all_means_groups==2))));
% Mean trial-rest
m_plot_trialrest =  tanh(mean(atanh(all_means(all_means_groups==3))));
CI_plot_trialrest = tanh(mean(atanh(all_CI95(all_means_groups==3))));

%
% Plot means 
hold all
mean_trials_plot    = scatter(m_plot_trial, CI_plot_trial ,40, 'filled', 'd',   'MarkerEdgeColor', 'k', 'MarkerFaceColor', colours(1, :) , 'LineWidth',1);
mean_rest_plot      = scatter(m_plot_rest, CI_plot_rest, 40, 'filled', 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor',  colours(2, :), 'LineWidth',1);
mean_trialrest_plot = scatter(m_plot_trialrest, CI_plot_trialrest, 40, 'filled', 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colours(3, :), 'LineWidth',1);

[~, objh] = legend(plot_mean_CI);
objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 5); %// set marker size as desired
%%
name_fig = get(gcf, 'name');



grid on
ylabel('CI95(RHOs) abs value');
xlabel('mean(RHOs)');

% plot max/min mean + CI95% (for eyeballing)
subplot(2,1,2);
[~, sort_ind] = sort(all_means(:));
% Size marker
sz = 2;
p(1) = scatter(all_means(sort_ind), all_means(sort_ind)+all_CI95(sort_ind), sz,  'filled');
hold all;
grid on
p(2) = scatter(all_means(sort_ind), all_means(sort_ind)-all_CI95(sort_ind), sz,  'filled');
p(3) = plot(all_means(sort_ind), all_means(sort_ind), ':');

for i = 1:2
    p(i).MarkerFaceAlpha = 0.1;
end 

legend({'mean(RHO)+CI95'; 'mean(RHOs)-CI95'; 'mean(RHOs)'});
set(legend, 'Location', 'southeast');
ylabel('CI95 max/min (RHOs)');
xlabel('mean(RHOs)');


%% Add suptitle
if strcmp(data_type, 'wholebrain')|| strcmp(data_type, 'empirical')
    voxelinfo = sprintf('whole brain - n_voxels[min=%i,max=%i,mean=%g]', voxelrange(1), voxelrange(2), round(mean_n_voxels));
elseif strcmp(data_type, 'ROI')
    voxelinfo = sprintf('ROI - n_voxels[min=%i,max=%i,mean=%g]', voxelrange(1), voxelrange(2), round(mean_n_voxels));
elseif strcmp(data_type, 'simulated')
    voxelinfo = sprintf('simulated - n_voxels[%i]', voxelrange(1));
end

suptitle(voxelinfo);

%% Save
%Print in command window
layout_indent(10); fprintf('Plotting of mean vs CI95 -- done'); layout_line_break(2);
figname = fullfile(targetdir, sprintf('fig-%s', get(gcf, 'name')));
if save_figures == 2
    disp(['Storing: ' figname]);
    save_fig(figname); % without TDT: save_fig, but then only .fig (check help for .eps and .png)
end
layout_line_break(1);

end