%% function mask_names = create_masks(cfg, nii_nam, nii_dat, subj)
%
% Function that creates masks for each model specified by
% generate_nifti_from_models. nii_dat and nii_nam are of run-size.
%
% Made by: Kai/Leonardo, 2019

% Updated to include phase-general masks independent of class - Ingmar
% 07-01-2021
% Separated from generate_nifti_from_models, combined with create_mask.m,
% added multiple other mask options - Ingmar 11-06-2020

function [cfg, mask_names] = create_masks(cfg, nii_nam, nii_dat, subj)
%% From cfg

preproc_dir = cfg.preproc_dir;
preproc_sub = fullfile(preproc_dir, sprintf('sub-%02i', subj));
phase       = cfg.saa.data_gen.model_names;

% class names
if isfield(cfg.saa, 'class_nodiff')
    class_str = cfg.saa.class_nodiff; 
else
    class_str = {'Aeh', 'Kon', 'Sum'}; % Sim Con/Rot Add
end 

%% Get description of data (should be the same for each image)

for d4_ind = 2:size(nii_nam, 4)
    assert(isequal(nii_nam(:, :, :, 1), nii_nam(:, :, :, d4_ind)), ...
        'generate_nifti_from_models:nii_nam_not_equal_across_time', ...
        'nii_nam_not_equal_across_time: the voxel description (nii_nam) differs (at least) between image 1 and %i', d4_ind);
end
nii_model_descr = nii_nam(:,:,:,1);

%% Generate masks

full_mask   = false(size(nii_dat(:,:,:,1))); % init
m_count     = 1; % init

for l_ind = 1:size(nii_dat, 3)
    
    % Initialize standard_mask
    std_mask            = false(size(nii_dat(:,:,:,1)));
    std_mask(:,:,l_ind) = ~cellfun(@isempty, nii_model_descr(:, :, l_ind)); % cellstr
    
    %% One mask for each full model (stored in the third dimension)
    
    description_mask    = sprintf('mask_fullmodel_%i', l_ind);
    fname_mask          = fullfile(preproc_sub, sprintf('sub-%02i_mask_fullmodel_%i.nii', subj, l_ind));
    mask_names{m_count} = fname_mask; m_count = m_count + 1;
    mask                = std_mask; % get mask for current model (currently: one mask per slice)
    
    % Assert that the mask does not contain nans
    assert(any(any(~any(isnan(full_mask)))), 'generate_nifti_from_models:mask_nans_present', ...
        'Potential mask %i (%s) contains %i NaNs for some reason. Aborting.', ...
        l_ind, description_mask, sum(isnan(full_mask)));
    
    % Only create a mask if the mask contains at least one true voxel
    if any(mask(:))
        create_mask(mask, fname_mask, description_mask);
        full_mask = full_mask | mask; % combine to full mask
    else
        warning('generate_nifti_from_models:skipping_empty_mask', ...
            'Data (%i/%i), fname: %s descr: %s, does not contain any true voxels. Skipping writing it.', ...
            l_ind, size(nii_dat, 3), fname_mask, description_mask);
        mask_names{l_ind} = '';
    end
            
    % Filter out any empty cells of model names, set to "empty"
    curr_nii_nam                = nii_nam(:,:,l_ind,:);
    empty_cells                 = cellfun(@isempty, nii_nam(:,:,l_ind,:));
    curr_nii_nam(empty_cells)   = {'empty'};
    nii_nam(:,:,l_ind,:)        = curr_nii_nam;
    
    %% Generate other masks, e.g. condition-specific, task trials, mixed effect voxels
    
    % Add padding
    zeros_padding   = zeros(size(full_mask(:, :, 1:end-1)));
    all_layers      = 1:size(nii_dat, 3);
    padding_dim     = all_layers(all_layers ~= l_ind);
    
    %% Masks for each phase independent of task class
    
    for ph_ind1 = 1:length(phase) % 3 phases
        
        warningv('create_mask:skip_cue', ...
            'Skipping cue-phase masking. Enable in create_masks.m if desired.');
        if strcmpi(phase(ph_ind1), 'modelcue'), continue; end
        
        % Get mask data
        m_descr_phase_only = sprintf('mask_%s_%i', phase{ph_ind1}, l_ind);
        fname_m_ph_only = fullfile(preproc_sub, sprintf('sub-%02i_%s.nii', ...
            subj, m_descr_phase_only));
        mask_names{m_count} = fname_m_ph_only; m_count = m_count + 1;
        
        % Assumes condition is stored at beginning of string
        ph_filter(:,:,l_ind) = contains(nii_nam(:,:,l_ind,1), ...
            phase{ph_ind1}, 'IgnoreCase', true); %#ok<*AGROW>
        ph_filter(:,:,padding_dim) = zeros_padding;
        create_mask(ph_filter, fname_m_ph_only, m_descr_phase_only);
    end
       
    %% Mask phase and task dependent masks
    
    % Check if phase-only masks should be created or full masks
    if ~isfield(cfg.saa, 'phase_masks_only') || ~cfg.saa.phase_masks_only
        
        % Check if all masks should really be created
        if ~isfield(cfg.saa, 'all_masks') || ~cfg.saa.all_masks
            continue; 
        else
            cfg.saa.all_masks = 1; 
        end
    else
        continue;
    end

    % Masks for each condition (block phases) 
    for cl_ind = 1:length(class_str) % 3 task classes
        
        m_descr_cl = sprintf('mask_%s_%i', class_str{cl_ind}, l_ind);
        fname_m_cl = fullfile(preproc_sub, sprintf('sub-%02i_%s.nii', subj, m_descr_cl));
        mask_names{m_count} = fname_m_cl; m_count = m_count + 1;
        
        % Assumes condition is stored at beginning of string
        t_filter(:, :, l_ind) = strncmpi(nii_nam(:,:, l_ind, 1), class_str{cl_ind}, 3); %#ok<*AGROW>
        t_filter(:, :, padding_dim) = zeros_padding;
        create_mask(t_filter, fname_m_cl, m_descr_cl);
        
        % Generate masks for task trials 3x3 (model phase)
        for ph_ind2 = 1:length(phase)
            
            % eg sub-02_mask_modelcue_Aeh_1
            m_descr_ph = sprintf('mask_%s_%s_%i', class_str{cl_ind}, phase{ph_ind2}, l_ind);
            fname_m_ph = fullfile(preproc_sub, sprintf('sub-%02i_%s.nii', subj, m_descr_ph));
            mask_names{m_count} = fname_m_ph; m_count = m_count + 1;
            
            % eg Sum_modelcue_fullmodel1 (variable str length)
            filter_str = sprintf('%s_%s_fullmodel%i', class_str{cl_ind}, phase{ph_ind2}, l_ind);
            t_ph_filter(:, :, l_ind) = strcmpi(nii_nam(:, :, l_ind, 1), filter_str);
            t_ph_filter(:, :, padding_dim) = zeros_padding;
            create_mask(t_ph_filter, fname_m_ph, m_descr_ph);
        end
    end
       
    %% Generate mask for mixed effects
    
    % Skip if no mixed voxels were added
    if ~isfield(cfg.saa.data_gen, 'n_mixedvoxels') || ...
            isempty(cfg.saa.data_gen.n_mixedvoxels) || ...
            all(cfg.saa.data_gen.n_mixedvoxels) == 0
        dispv(1, '    No mixed voxels, mask skipped');
    else        
        m_descr_mix = sprintf('mask_mixed_fx_%i', l_ind);
        fname_mix = fullfile(preproc_sub, sprintf('sub-%02i_%s.nii', subj, m_descr_mix));
        mask_names{m_count} = fname_mix; m_count = m_count + 1;
        
        mix_fx_filter(:, :, l_ind) = strncmpi('mixed', nii_nam(:, :, l_ind, 1), 5);
        mix_fx_filter(:, :, padding_dim) = zeros_padding;
        create_mask(mix_fx_filter, fname_mix, m_descr_mix);
    end
    
    % Print
    dispv(1, '\n    Masks for model %0i done (one layer)...\n', l_ind);
end

dispv(1, '    Completed mask generation for %0i models\n', l_ind);
mask_names = mask_names(~cellfun(@isempty, mask_names)); % remove empty masks

%% Assert that full mask covers all data

[dim1, dim2, ~, ~]      = size(nii_dat);
[maskdim1, maskdim2, ~] = size(full_mask);
assert(maskdim1 && maskdim2 == dim1 && dim2, 'generate_nifti_from_models:no_full_data_masking', ...
    'Full mask does not fully cover simulated voxel data. Aborting.');

%% Saving full mask

description_mask    = 'mask for all voxels of any simulation';
fname_mask          = fullfile(preproc_sub, 'mask.nii');
create_mask(full_mask, fname_mask, description_mask);

% Save metadata (masknames, names of voxels)
masknames_fname = fullfile(preproc_sub, 'mask_names.mat');
disp(['    Saving names of masks for different models to ' masknames_fname])
info  = [];
info.mask_names         = 'Name of maskfiles for different simulations';
info.nii_model_descr    = 'Names of simulation in each voxel (across time, first vol in nii_nam)';
save(masknames_fname, 'mask_names', 'nii_model_descr', 'info');

% Save description of each voxel
info.nii_nam = 'Names of simulation in each voxel for each timepoint (nii_nam)';
niinam_fname = fullfile(preproc_sub, 'nii_nam.mat');
disp(['    Saving  description of each voxel (nii_nam and nii_model_descr) to ' niinam_fname])
save(niinam_fname, 'nii_nam', 'nii_model_descr');

% Storing to cfg
cfg.mask_names      = mask_names;
cfg.nii_model_descr = nii_model_descr;

end



%% Subfunction: create_mask(mask, fname, description)
% 
% Function that creates a nifti image
% mask: mask as 3d matrix, containing false and trues, with dimensions of
%   data
% fname: name of the output
% description (optional): info about each voxel, it can be either a string
% or struct with two fields called 'name' and 'voxels'

function create_mask(mask, fname, description)

clear ni img hdr

% If file name does not have .nii extension, add it
if  ~contains(fname, '.nii')
    fname = [fname '.nii'];
end

% Create image
ni = nifti;

% Make header info
ni.dat       = file_array(fname, size(mask), [0 spm_platform('bigend')], 0, 1, 0); % SPM
ni.mat       = eye(4);
ni.descrip   = description; % CHANGE THIS?
ni.dat.dtype = 'FLOAT32';

% Compress name for printing info
if length(fname) > 80
    [p, n, x] = fileparts(fname);
    dispv(2, '    Writing: ...%s%s%s%s', p(end-60:end),filesep, n, x);
else
    dispv(2, '    Writing: %s', fname);
end

create(ni); % SPM func

% write image
ni.dat(:,:,:) = mask(:,:,:);

end
