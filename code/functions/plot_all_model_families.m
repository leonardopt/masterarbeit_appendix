% Function to plot and store all model plots, convolved with the hrf

% requires hrf file in the folder
% figure output given by user


function plot_all_model_families(sim_name, targetdir, dt, response)

set(0, 'DefaultFigureVisible', 'off')

if ~exist('dt', 'var')
    dt = 0.0606;
end

if ~exist('response', 'var')
    response = [1 1 1];
end


if ~exist('sim_name', 'var')
    sim_name     = {'sim-boxcar',...
        'sim-sticknoaftereffect_onestickrest', ...
        'sim-sticknoaftereffect_allstickrest', ...
        'sim-sticknoaftereffect_boxrest', ...
        'sim-stickfadeoff', ...
        'sim-sticksuppression',...
        'sim-stickshortsuppression',...
        'sim-piecewisesuppression'....
        'sim-piecewisefadeoff',...
        };
end



for h_ind = 1:length(sim_name)
    
    [modelcue, modeltrial, modelrest] = models_from_hypothesis(sim_name{h_ind}, dt, response, 0);
    
    simfolder = fullfile(targetdir,'model_plots',  sim_name{h_ind} );
    [a,b ] = mkdir(simfolder);
    
    modelname = sim_name{h_ind}(5:end);
    cueplot = plot_models_and_hrf( modelcue, modelname) ;
    cuefigname = fullfile(simfolder ,'model_cue');
    save_fig(cuefigname);
    
    trialplot =  plot_models_and_hrf( modeltrial,modelname) ;
    trialfigname = fullfile(simfolder ,'model_trial');
    save_fig(trialfigname);
    
    restplot =  plot_models_and_hrf( modelrest,modelname) ;
    restfigname = fullfile(simfolder ,'model_rest');
    save_fig(restfigname);
    
    
end

set(0, 'DefaultFigureVisible', 'on')

end