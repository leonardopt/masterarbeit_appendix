% Meta caller script to run the runfile multiple times, adjusting only some
% essential parameters. This specific script calls all specified
% simulations with different parameters.
%
% NB1: Most parameters are not looped here but passed as [] which the
% underlying data generation functions will take one at a time in multiple
% for loops. Therefore, it does not make sense to pass multiple {[] [] []}.
% It is possible to loop multiple hypotheses (so aftereffect models, not to
% be confused with the timeseries models per layer).
%
% NB2: Moreover, as these functions are actually optimized to create only
% one dataset per "subject", not multiple as batch, it should be passed
% carefully and perhaps checked along the way, to see if indeed all
% parameter sets survive (and not just the last one is taken).
%
% NB3: Lastly, it
% is currently not possible/feasible to pass different voxel dimensions per
% layer, so e.g. changing the number of mixed voxels. Perhaps just set the
% noise levels of those to 0 or maybe program them to be NaNs after
% resizing data in create_nifti_multiple_models.m
%
% Leo 10-11

function cfg = meta_runfile_sim_control_mixedeffectvoxels5

scriptdir = pwd;

%% Do struct

do.modelgeneration          = 0;
do.firstlevel               = 0;
do.decoding                 = 0;
do.decoding_grouplevel      = 0; % group level decoding 
do.parfor                   = 1; % turn off for bugfixing!
do.all_subs                 = 1;
do.similarityanalysis       = 0;
do.similarityanalysis_group = 1; % group level similarity
% Store in cfg struct
cfg.do = do;

%%  Display output figures
% Comment out if you want to show figures a long as they're created
set(0, 'DefaultFigureVisible', 'off')
% Set some defaults for plotting
set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(0, 'DefaultAxesTickLabelInterpreter', 'none')

%% Settings
% Set paths
cfg.all_paths   = set_all_paths(); % global_all_paths_private
cfg.subs_todo   = get_subs_list() 
cfg.data_type   = 'simulated';
% Folder names
cfg.level1      = 'level1-FIR1';
cfg.post_level1 = 'decoding-timeresolved_multiclass-cv';
cfg.derivative_dir  = fullfile(cfg.all_paths.data_base_dir, 'derivatives_test_sim_control_mixedeffectvoxels5');
cfg.saa.data_gen.batch_version = 'Leo-27-11-20'; % update if needed
% Cross/xclass validation
cfg.analysis_type = 'standard_cv';
% Condition/phase/block "models"
cfg.saa.data_gen.models.weightformat    = 'matrixcell'; % might be redundant
cfg.saa.data_gen.dt                     = 2/33;
cfg.saa.data_gen.model_names            = {'modelcue', 'modeltrial', 'modelrest'};
% Kongruent, Summe, Aehnlich
cfg.saa.tasks = {'aeh-s', 'aeh-c', 'kon-s', 'kon-c', 'sum-s', 'sum-c'}
% Convolution
cfg.saa.data_gen.hrf_conv = 'hrf_conv';
% Output
cfg.results.output = {'accuracy', 'confusion_matrix'};
% Other settings
cfg.saa.all_masks                   = 0; % mask every voxel type, turn off if not needed (only fullmodels)
cfg.saa.data_gen.plotting           = 0; % plot stick model and post-hrf timeseries
cfg.flags.overwrite_level1_SPM      = 0;
cfg.flags.data_gen_overwrite_niftis = 0;
cfg.parfor_workers                  = 24; 

% SPM automatic masking 
% 1 if you want to use SPM masking, 0 if you want to deactivate it
cfg.saa.spmautomaticmasking = 1;



%% ---- Parameters per 3rd dim (voxel layer) ---- %%

% Initialise layers 
layer_descr         = {'weights 1', ...
                      'weights 2', ...
                      'weights 3', ...
                       'weights 4', ...
                       'weights 5' ... 
                            };
layer_number        = [1:length(layer_descr)];

% Voxels (currently not possible to vary these per layer, see NB3)
task_voxels         = 9;
n_mixedvoxels       = [5 ...
                       5 ... 
                       5 ... 
                       5 ...
                       5 ...
                         ];
total_n_voxels      = task_voxels + n_mixedvoxels(1); 

% Weights (normalized below)
taskweights         = {  [0.4 1 0.7] ...
                         [0.4 1 0.7] ...
                         [1.2 2 1.6] ...
                         [1.2 2 1.6] ...
                         [1.5 2.5 2] ...
                          }; % Aeh Kon Sum
phaseweights        = {[0.1 0.6 0.2] ...
                       [0.1 0.6 0.2] ...
                       [0.1 0.6 0.2] ...
                       [0.1 0.6 0.2] ...
                       [0.1 0.6 0.2] ...
                        }; % Cue Trial Rest
weight_mixed_voxels = [0.5 ...
                       1 ... 
                       0.5 ... 
                       1 ...
                       1 ...
                        ];
priorweights        = {1.5 ...
                       1.5 ...
                       1.5 ...
                       1.5 ...
                       1.5 ...
                        }; % baseline weight
                   
% Weight normalization
%taskweights     = cellfun(@(x) x/sum(x), taskweights, 'UniformOutput', false);
%phaseweights    = cellfun(@(x) x/sum(x), phaseweights, 'UniformOutput', false);
%normalise_weight_mat = 1; % normalise weight matrix, so that each row sums up to 1

% Noise
perc_voxelnoise   = [0.01 ...
                     0.01 ...
                     0.01 ...
                     0.01 ...
                     0.01 ...
                     ];
noise_after_conv    = [0.01 ...
                       0.01 ...
                       0.01 ...
                       0.01 ...
                       0.01 ...
                        ];
perc_additionalnoise =  [0.1 ...
                         0.1 ...
                         0.1 ...
                         0.1 ...
                         0.1 ...
                         ];

% Check dimensions match
parameters = [length(n_mixedvoxels) ...
              length(taskweights) ...
              length(phaseweights) ...
              length(weight_mixed_voxels) ...
              length(priorweights) ...
              length(perc_voxelnoise) ...
              length(noise_after_conv) ...
              length(perc_additionalnoise) ...
              ];
          
if sum(parameters) ~= length(layer_number) * length(parameters)
    error('Parameters given of inconsisten number. Please check the layers!')
end 
      
      
% Condition/phase/block "model"
response            = [1 10 10]; % Cannot vary per layer! Use single set of responses per batch
restboxcar          = 0;
trialboxcar         = 0;

% Initialise table 
T = table();
T.layer_number      = layer_number';
T.layer_descr       = layer_descr';
T.task_voxels       = repmat(task_voxels, max(T.layer_number), 1);
T.n_mixedvoxels     = n_mixedvoxels';
T.total_n_voxels    = repmat(total_n_voxels, max(T.layer_number), 1);
T.model_names         = repmat(cfg.saa.data_gen.model_names, max(T.layer_number), 1);
T.taskweights         = cell2mat(taskweights');
T.phaseweights        = cell2mat(phaseweights');
T.weight_mixed_voxels = weight_mixed_voxels';
T.priorweights        = cell2mat(priorweights');
T.perc_voxelnoise = perc_voxelnoise';
T.noise_after_conv  = noise_after_conv';
T.perc_additionalnoise   = perc_additionalnoise';
T = addvars(T, repmat(response, max(T.layer_number), 1), 'After', 'model_names', 'NewVariableNames', {'response'});
T = addvars(T, repmat(cfg.saa.data_gen.dt,  max(T.layer_number), 1), 'After', 'response', 'NewVariableNames', 'dt');


%% Assign parameters to structure fields

% Parameters
cfg.saa.data_gen.restboxcar             = restboxcar;
cfg.saa.data_gen.trialboxcar            = trialboxcar;
cfg.saa.data_gen.total_n_voxels         = total_n_voxels;
cfg.saa.data_gen.task_voxels            = task_voxels;
cfg.saa.data_gen.n_mixedvoxels          = n_mixedvoxels;
cfg.saa.data_gen.weight_mixed_voxels    = weight_mixed_voxels;
cfg.saa.data_gen.taskweights            = taskweights;
cfg.saa.data_gen.phaseweights           = phaseweights;
cfg.saa.data_gen.priorweights           = priorweights;
cfg.saa.data_gen.perc_voxelnoise        = perc_voxelnoise;
cfg.saa.data_gen.noise_after_conv       = noise_after_conv;
cfg.saa.data_gen.perc_additionalnoise   = perc_additionalnoise;
cfg.saa.data_gen.response               = response;
cfg.saa.data_gen.layer_descr            = layer_descr;
cfg.saa.data_gen.layer_number           = layer_number;

% Voxel and models setup
for p_ind = 1:length(taskweights)
    % ---------------------------------------------------------------%
    pref_weight{p_ind} = define_voxel_weights_matrix(priorweights{p_ind}, phaseweights{p_ind}, taskweights{p_ind});
    % ---------------------------------------------------------------%
    

end
cfg.saa.data_gen.pref_weight = pref_weight; % Place


%% ---- Hypotheses/models ---- %%

sim_name     = {...
                 'sim-boxcar',...
                  'sim-sticknoaftereffect_allstickrest', ...
                'sim-sticknoaftereffect_randstickrest' ...
                 'sim-sticknoaftereffect_boxrest', ...
                'sim-stickfadeoff', ...
                'sim-sticksuppression',...
                'sim-stickshortsuppression',...
                'sim-piecewisesuppression'....
               'sim-piecewisefadeoff',...
                  'sim-sticknoaftereffect_randstickrest_higherfirststick', ...
                };

% % Create table with all simulations
% S = table();
% for s = 1:length(sim_name)
%     curr_T = T;
%     sim = cellstr(repmat(sim_name{s}(5:end), max(T.layer_number), 1));
%     curr_T = addvars(curr_T, sim, 'Before', 'layer_number');
%     S = vertcat(S, curr_T);
% end

%% Save tables
[a,b ] = mkdir(cfg.derivative_dir);
writetable(T, fullfile(cfg.derivative_dir, 'sim_shared_parameters.csv'))
           
%% Run full pipeline for each hypothesis 

n_hypotheses = length(sim_name);
cfg.saa.data_gen.aftereffect = '';
for h_ind = 1:n_hypotheses
    
     cfg.sim_name = sim_name{h_ind};
    %cfg.saa.data_gen.aftereffect = aftereffects{h_ind};
        
    % Only for piecewise defined functions: reduce number of voxels to 3
    % (only task related)
    if contains(sim_name{h_ind}, 'piecewise')
        cfg.saa.data_gen.ninetothreevoxels = 1;
    else
        cfg.saa.data_gen.ninetothreevoxels = 0;
    end
    
    % Phase models ------------------------------------------------------%
    [modelcue, modeltrial, modelrest] = models_from_hypothesis(sim_name{h_ind}, cfg.saa.data_gen.dt, response, 0);
    % -------------------------------------------------------------------%
    cfg.saa.data_gen.all_models = {modelcue modeltrial modelrest};
    
    % Recalculate noise parameters as a percentage of model heights
    maxval_before_conv = max(cellfun(@max, cfg.saa.data_gen.all_models));
    % Noise before convolution
    noise_before_conv = maxval_before_conv * perc_voxelnoise;
%     noise_before_conv =  perc_voxelnoise;

    % Additional noise
    additionalnoise = maxval_before_conv * perc_additionalnoise;
%     additionalnoise =  perc_additionalnoise;

    cfg.saa.data_gen.noise_before_conv = noise_before_conv;
    cfg.saa.data_gen.additionalnoise   = additionalnoise;
    
    curr_T = T;
    curr_T.noise_before_conv = noise_before_conv';
    curr_T.additionalnoise   = additionalnoise';
    % Add name of current simulation
    function_set = cellstr(repmat(sim_name{h_ind}(5:end), max(T.layer_number), 1));
    curr_T = addvars(curr_T, function_set, 'Before', 'layer_number');
    % Save table 
    writetable(curr_T, fullfile(cfg.derivative_dir, [sim_name{h_ind} '_parameters.csv']))

    % --------------------------------- Run ---------------------------- %
    run_timeresolved_analyses(cfg, do);
    % -------------------------------------------------------------------%
    
    %% Reduce SPM mat to free space
    disp('Reduce SPM mat size to save space...')
    reduce_SPM_mats_from_dirs(cfg)

    
    cd(scriptdir)
end % aftereffect hypotheses

cd(cfg.derivative_dir)
% %% Store all models 
% for s = 1:length(sim_name)
%     simulation_dir = fullfile(cfg.derivative_dir, sim_name{s})
%     create_all_models_plots(simulation_dir, 1) 
% end 


end
