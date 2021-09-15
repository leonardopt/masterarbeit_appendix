% [cfg, layer_descr] = get_masks_to_decode(cfg, sub_sourcedir)
%
% Split up from main body of afx_tdt_subject_timeresolved_decoding.m
%
% Patched an occurence where mask_names was not copied to the level1 folder
% due to analysis controls, now also looks at original data
% folder.
%
% Ingmar 17-03-21

function [cfg, layer_descr] = get_masks_to_decode(cfg, sub_sourcedir)
%% Set input and output to cfg: masks
% goal: set cfg.files.masks
% read files from spm_mat, spm_reduced mat, or as fallback, from betas (or betas_cached)

cfg.analysis = cfg.data_type;

switch cfg.data_type
    case 'simulated'
        try
            mask_names_mat = load(fullfile(sub_sourcedir, 'mask_names.mat')); % get mask names
        catch e1
            disp(['    Error loading: ' e1.message]);
            try
                preproc_subdir = fullfile(cfg.preproc_dir, sprintf('sub-%02i', cfg.subj), 'mask_names.mat');
                mask_names_mat = load(preproc_subdir); % get mask names
            catch e
                disp(['    Error: ' e.message]);
                error('Could not load mask_names.mat. Check mask creation, input, and file locations.');
            end
        end
        cfg.files.mask  = mask_names_mat.mask_names;
        cfg.analysis    = 'ROI'; % one 'ROI' analysis per mask
        display(cfg.files.mask)
        
        cfg.save_figures            = 2;
        cfg.plot_design             = 0; %1: show and save, 2: save only, 0: dont show or save
        cfg.plot_selected_voxels    = 0;  % check all masks are unique
        
    case 'ROI'
        
        if ~isfield(cfg, 'files') || ~isfield(cfg.files, 'mask') || isempty(cfg.files.mask)
            
            % Try this, because cfg.files.mask is otherwise dependent on
            % sourcedir (cfg.level1_dir) even though masks could be elsewhere
            errstring = ['No mask name input found for ROI analysis, aborting analysis. ', ...
                'Specify derivative_dir subfolder for native-space ROI masks in cfg.ROI_folder, ', ...
                'and the names in cfg.ROI_masks. Assumes 1 folder per participant, ', ...
                '"subjectspace_ROIs", with >=1 mask nifti. ',...
                'Alternatively, pass full filenames directly to cfg.files.mask.'];
            assert(isfield(cfg, 'ROI_masks') && isfield(cfg, 'ROI_folder'), ...
                errstring);
            mask_loc = fullfile(cfg.derivative_dir, cfg.ROI_folder, ...
                sprintf('sub-%02i', cfg.subj), 'subjectspace_ROIs');
            dispv(1, '    Selecting masks from: %s', mask_loc);
            cfg.files.mask = fullfile(mask_loc, cfg.ROI_masks); % init
            assert(all(isfile(cfg.files.mask)), 'ROI mask files are missing, please check!');
        else
            dispv(1, '    Current ROI mask passed directly via cfg.files.mask:'); disp(cfg.files.mask);
        end
        
        cfg.save_figures    = 2;
        cfg.plot_design     = 0;
        
    case {'wholebrain', 'searchlight'}
        
        % nothing to do, will be handled by TDT
        cfg.files.mask      = fullfile(sub_sourcedir, 'mask.nii');
        cfg.analysis        = cfg.data_type;
        cfg.save_figures    = 0;
    otherwise
        error('Unknown datatype for masks')
end

% get layer descriptions if simulation
layer_descr = [];
try layer_descr = cfg.saa.data_gen.layer_descr; end
        
end
