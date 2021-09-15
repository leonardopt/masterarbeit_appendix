% [regressor_names, beta_fnames] = load_and_validate_betas_for_decoding(sub_sourcedir)
%
% Split up from main body of afx_tdt_subject_timeresolved_decoding.m to
% re-use for searchlight analysis.
%
% Ingmar 17-03-21

function [regressor_names, beta_fnames] = load_and_validate_betas_for_decoding(sub_sourcedir)

% Get betas in folder
dispv(1, '1. Getting betas in folder...');

% Use at least for verification of regressor information below
beta_fnames = cellstr(spm_select('FPList', sub_sourcedir, '^beta_[0-9]*.nii'));
assert(~isempty(beta_fnames), 'No beta images found in %s, aborting', sub_sourcedir);

try % create regressor_names.mat
    if isfile(fullfile(sub_sourcedir, 'regressor_names.mat')) 
        % Load data from regressor_names.mat
        disp(['  Getting regressors from ' fullfile(sub_sourcedir, 'regressor_names.mat')]);
        load(fullfile(sub_sourcedir, 'regressor_names.mat'), 'regressor_names');
        % SPM.mat exist: create regressor names from SPM.mat
    elseif isfile(fullfile(sub_sourcedir, 'SPM.mat')) 
        disp(['  Getting regressors from SPM.mat and creating ' fullfile(sub_sourcedir, 'regressor_names.mat')]);
        regressor_names = design_from_spm(sub_sourcedir); % will also create regressor_names.mat
    end
catch
    warning('Getting regressor names from regressor_names.mat or SPM.mat was not successful')
    regressor_names = '';
end

if isempty(regressor_names)
    warning(['Getting regressor names from regressor_names.mat or SPM.mat was not successful, ', ...
        'trying fallback option to generate regressor_names.mat from images'])
    if ~isempty(beta_fnames) % Betas exist: create regressor_names.mat from betas
        dispv(1, '  Getting regressors from beta files...');
        hdrs = spm_vol(beta_fnames);
        regressors = cell(size(hdrs));
        
        % Extract name from description, example: 'spm_spm:beta (2167) - Sn(8) AehCRest 32*bf(1)' => Sn(8) AehCRest 32
        for h_ind = 1:numel(regressors)
            parts = strsplit(hdrs{h_ind}.descrip, {' - '});  % split at ' - '
            regressors{h_ind} = parts{2}; %only keep part 2
        end
        
        % Create regressor_names.mat
        regressor_names = design_from_spm(regressors,'regressors');
        disp([' Getting regressors from beta files successful, writing ' fullfile(sub_sourcedir, 'regressor_names.mat')])
        save(fullfile(sub_sourcedir, 'regressor_names.mat'), 'regressor_names');
    end
end

if isempty(regressor_names)
    error('Cannot identify input files. No regressor_names.mat, SPM.mat, or beta files detected in %s, or extraction from these sources did not work', source_dir)
end

%% Loaded betas checks

assert(isfile(fullfile(sub_sourcedir, 'regressor_names.mat')), ...
    'error regressor_names.mat does not exist, should have exist already or been created. Please check.')

% Check that regressor information is up to date & check number of betas
assert(size(beta_fnames, 1) == size(regressor_names, 2), ...
    'Number of betas differs from number of regressors, maybe data regressor names is corrupt. Please check.');

% Check first, last, and few random headers
selected_idx = sort([1, randperm(length(beta_fnames)-2, 4)+1, length(beta_fnames)]);
disp(['2. Data consistency check: Checking regressors vs headers of betas for regressors: ' num2str(selected_idx)]);
selected_hdrs = spm_vol(beta_fnames(selected_idx, :)); % SPM check

for h_ind = 1:numel(selected_hdrs)
    assert(contains(selected_hdrs{h_ind}.descrip, regressor_names{3, selected_idx(h_ind)}), ...
        ['Data consistency check failed for regressor %i: no match between', ...
        ' regressor_names{3, %i}="%s" and .descrip in file (maybe field does not exist?).\nFilename: %s'], ...
        selected_idx(h_ind), selected_idx(h_ind), regressor_names{3, selected_idx(h_ind)}, selected_hdrs{h_ind}.fname)
end

dispv(1, '  Data consistency check: successful');


end % func