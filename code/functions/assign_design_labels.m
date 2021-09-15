% function all_cfg = assign_design_labels(cfg, regressor_names, sourcedir)
%
% Describe data (set labels)
% 'standard_cv' 'non_xclass_difficulty_generalization' 'xclass_difficulty_control'
% 'crosstask_difficulty_generalization' 'difficulty_generalization'
%
% Dependent on afx_tdt_subject_timeresolved_decoding.m
%
% Kai and Ingmar, 06/2020

% Added paths - Ingmar 08/01/2021

function all_cfgs = assign_design_labels(cfg, regressor_names, sourcedir)

dispv(1, '3. Design construction... \n  Labeling data according to analysis type "%s".', cfg.analysis_type);

if ~isfield(cfg, 'software') || isempty(cfg.software) || strcmp(cfg.software, 'none')
    cfg.software = 'SPM12'; % set
    try
        global_all_paths();
        addpath(genpath(all_paths.spm12)); addpath(genpath(all_paths.tdt));
    catch
        dispv(1, '    SPM12 was not found on path, it might stop design construction. Add path manually or check global_all_paths if it does.');
    end
end

all_cfgs = {}; % init

% Typo
if strcmpi(cfg.analysis_type, 'standard'), cfg.analysis_type = 'standard_cv'; end

switch cfg.analysis_type
    
    case 'pairwise_cv' 
        % ----------------------------------------------------------------------%
        % Used for searchlight decoding (the default pipeline uses a different
        % pair-wise method of standard_cv below)
        labelnames  = {'Aeh*', 'Sum*'} % repeat for other label combinations
        labels      = [1,     2]
        all_cfgs{1} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir);
        % ----------------------------------------------------------------------%
        labelnames  = {'Aeh*', 'Kon*'}
        all_cfgs{2} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir);
        % ----------------------------------------------------------------------%
        labelnames  = {'Kon*','Sum*'}
        all_cfgs{3} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir);
        % ----------------------------------------------------------------------%
        
    case 'standard_cv'
        % ----------------------------------------------------------------------%
        % Simplest comparison (unbalanced), difficulty levels disregarded
        labelnames  = {'Aeh*', 'Sum*', 'Kon*'}  %#ok<*NOPRT>
        labels      = [1,       2,      3]
        all_cfgs{1} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir);
        % ----------------------------------------------------------------------%
        
    case 'standard_difficulty_generalization'
        % ----------------------------------------------------------------------%
        % Simplest comparison of difficulty levels, task conditions
        % disregarded
        labelnames      = {'AehS*', 'AehC*', 'SumS*', 'SumC*', 'KonS*', 'KonC*'}
        labels          = [1,       2,      1,      2,      1,      2]
        all_cfgs{1} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir);
        % ----------------------------------------------------------------------%
        
    case 'xclass_difficulty_generalization'
        % ----------------------------------------------------------------------%
        % Cross-class decoding between task sets for difficulty
        % generalization across sets
        %    - original labelling intact (unbalanced)
        %    - 'one-way' since class information is not directly contrasted
        labelnames      = {'AehS*', 'AehC*', 'SumS*', 'SumC*', 'KonS*', 'KonC*'}
        labels          = [1,       1,      2,      2,      3,      3]
        xclass_labels   = [1,       2,      1,      2,      1,      2]
        all_cfgs{1} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir, xclass_labels);
        % ----------------------------------------------------------------------%
        warning('Currently, xclass is directional from simple to complex!');
        
    case 'xclass_difficulty_control'
        % ----------------------------------------------------------------------%
        % Control for task difficulty effects using one task set simple & complex vs. another,
        % but differentiating with a cross-class train and test switch.
        %    - train on task 1 regardless of difficulty, test on task 2
        %    - cross design to separate effects from 1 -> 2 vs. 2 -> 1
        %    - thus, directional training effects of difficulty are equalized/controlled
        cfg.files.twoway = 1; % train and test on all possible pairs
        labelnames      = {'AehS*', 'AehC*', 'SumS*', 'SumC*'} % repeat for other label combinations
        labels          = [1,       1,      2,      2]
        xclass_labels   = [1,       2,      2,      1]
        all_cfgs{1} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir,xclass_labels);
        % ----------------------------------------------------------------------%
        labelnames  = {'AehS*', 'AehC*', 'KonS*', 'KonC*'}
        all_cfgs{2} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir,xclass_labels);
        % ----------------------------------------------------------------------%
        labelnames  = {'KonS*', 'KonC*', 'SumS*', 'SumC*'}
        all_cfgs{3} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir,xclass_labels);
        % ----------------------------------------------------------------------%
        
    case 'crosstask_difficulty_generalization'
        % ----------------------------------------------------------------------%
        % Generalizes task difficulty patterns across tasks without
        % removing (task-specific) intraclass information ('real' xclass decoding)
        %    - train on simple, test on complex across tasks
        %    - then, separate task information using xclass labels (design)
        %    - non-directional since all pairs are trained and tested
        cfg.files.twoway = 1; % train and test on all possible pairs
        labelnames      = {'AehS*', 'AehC*', 'SumS*', 'SumC*'} % repeat for other label combinations
        labels          = [1,       2,      1,      2]
        xclass_labels   = [1,       1,      2,      2]
        all_cfgs{1} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir,xclass_labels);
        % ----------------------------------------------------------------------%
        labelnames  = {'AehS*', 'AehC*', 'KonS*', 'KonC*'}
        all_cfgs{2} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir,xclass_labels);
        % ----------------------------------------------------------------------%
        labelnames  = {'KonS*', 'KonC*', 'SumS*', 'SumC*'}
        all_cfgs{3} = decoding_describe_data(cfg,labelnames,labels,regressor_names,sourcedir,xclass_labels);
        % ----------------------------------------------------------------------%
        
    otherwise
        error('Unknown analysis type %s', cfg.analysis_type)
end
pause(1);