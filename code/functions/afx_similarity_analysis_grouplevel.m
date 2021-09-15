%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform similarity analysis at the group level
% 1) Load data
% 2) Average complex and simple trials
% 3) Plot mean correlation between trials and their CI95; plot CI92 vs mean
% 4) Plot time-resolved correlation
% 5) Save data (cfg file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function afx_similarity_analysis_grouplevel(run_cfg)
%% Process control
% removed do struct since that causes parfor issues

do_plotting = 1;
do_saving   = 1;

%% Get info

subs_todo = run_cfg.subs_todo;
if strcmp(run_cfg.data_type , 'ROI')
    sourcefolder = run_cfg.corr_dir_ROI;
else
    sourcefolder = run_cfg.corr_dir;
end

targetfolder = fullfile(sourcefolder, 'group');
if ~isfolder(targetfolder)
    mkdir(targetfolder)
else
    rmdir(targetfolder, 's')
    mkdir(targetfolder)
end

%% 1) Load data for all subjects

dispv(1, 'I. Starting data retrieval for all given subjects (n=%0i)', length(subs_todo));

all_similarity_cfg = cell(1, length(subs_todo));
layer_cfg = struct;

for s_ind = 1:length(subs_todo)
    subjdir = fullfile(sourcefolder, sprintf('sub-%02i', subs_todo(s_ind)));
    similarity_cfg = load(fullfile(subjdir, 'similarity_cfg.mat'));
    similarity_cfg = similarity_cfg.similarity_cfg;
    all_similarity_cfg{s_ind} = similarity_cfg;
end

%% SIMULATED

if strcmp(run_cfg.data_type, 'simulated')
    %% Reorganise the data
    
    for layer = run_cfg.saa.data_gen.layer_number
        for s_ind = 1:length(subs_todo)
            similarity_cfg = all_similarity_cfg{s_ind};
            layer_cfg(layer).all_subj_cfg{s_ind} = similarity_cfg(layer);
        end
    end
    
    %% Analyise model by model
    
    for layer = run_cfg.saa.data_gen.layer_number
        
        sprintf('-- ANALYSE LAYER %i --', layer);
        
        %% Get data from current model
        all_cfg_corr = layer_cfg(layer).all_subj_cfg;
        
        %% Average coplex and simple trials
        cfg_sim = afx_merge_complex_simple(all_cfg_corr);
        % Find and store some information about the total number of voxels
        disp('Average complex and simple - completed')
        
        %% Plot means and 95% confidence intervals
        % Update simualtion name and target directory
        if do_plotting
            cfg_sim.general_info.save_figures = 2; % 2 if you want to save them
            layerfolder = fullfile(targetfolder, sprintf('layer-%i', layer));
            if ~isfolder(layerfolder), mkdir(layerfolder); end
            cfg_sim.general_info.target_directory = layerfolder;
            cfg_sim = afx_plot_means_and_CIs(cfg_sim);
            disp('Plot means and CIs - completed')
            close all
            
            %% Plot time-resolved correlation
            % check "analysis_defaults" to invert even vs odd order
            % check "analysis_defaults" to plo single beta/average beta (default)
            % cfg_sim = afx_set_timeresolved_plotting_order(cfg_sim);
            
            % Plot all timecourses
            cfg_sim = afx_plot_all_timecourses(cfg_sim);
            % Plot averaged timecourse
            cfg_sim = afx_plot_averaged_timecourses(cfg_sim);
            disp('Plot time-resolved correlation - completed')
        end
    end
    
    % Saving results in
    if do_saving
        sprintf('Saving results in: %s', fullfile(targetfolder, 'cfg_sim'))
        save(fullfile(targetfolder, 'cfg_sim'), 'cfg_sim');
    end

else
    %% EMPIRICAL
    
    %% Average coplex and simple trials
    cfg_sim = afx_merge_complex_simple(all_similarity_cfg);
    % Find and store some information about the total number of voxels
    disp('Average complex and simple - completed')
    
    %% Plot means and 95% confidence intervals
    % Update simualtion name and target directory
    
    if do_plotting
        cfg_sim.general_info.save_figures = 2; % 2 if you want to save them
        cfg_sim.general_info.target_directory = targetfolder;
        cfg_sim = afx_plot_means_and_CIs(cfg_sim);
        disp('Plot means and CIs - completed')
        close all
        
        %% Plot time-resolved correlation
        % check "analysis_defaults" to invert even vs odd order
        % check "analysis_defaults" to plo single beta/average beta (default)
        % cfg_sim = afx_set_timeresolved_plotting_order(cfg_sim);
        
        % Plot all timecourses
        cfg_sim = afx_plot_all_timecourses(cfg_sim);
        % Plot averaged timecourse
        cfg_sim = afx_plot_averaged_timecourses(cfg_sim);
        disp('Plot time-resolved correlation - completed')
    end
    
    if do_saving
        % Saving results in
        sprintf('Saving results in: %s', fullfile(targetfolder, 'cfg_sim'))
        save(fullfile(targetfolder, 'cfg_sim'), 'cfg_sim');
    end
end