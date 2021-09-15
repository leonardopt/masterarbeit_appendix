%% function similarity_analysis_results = afx_compute_subject_means(cfg)
%
% Performs Fisher-Z transformed correlations on odd versus even-run beta
% images from a participant or a simulated voxel space.
%
% IN
%   cfg: dirs for level 1 and 2; cfg.subj = subject number
%       data_type = 'simulation'/'wholebrain'/'ROI'
%       simulation name if available
%       flags for optional control settings
%       optional design settings (e.g. make_pairwise)
%
% NEEDS: (masked) data, full list of timeresolved betas in sourcedir
%
% OUT
%   similarity_analysis_results: including all files and design matrices, 
%   correlation results, labels
%
% Kai and Leonardo, 2019

function similarity_analysis_results = afx_compute_subject_means(cfg, layer)

if ~exist('layer', 'var')
    layer = 'no layer';
end

%% Source and target folder
dispv(1, 'Starting afx similarity for cfg.subj %i', cfg.subj)

sourcedir = cfg.level1_dir; % e.g. /analysis/share/corinna_kai/afx/derivatives/level1-FIR1
% (/func)/sub-99 will be added below
if strcmp(cfg.data_type,'ROI')
    targetdir = cfg.corr_dir_ROI;
else
    if isfield(cfg, 'corr_dir')
        targetdir = cfg.corr_dir;
    else
        targetdir = fullfile(cfg.derivative_dir, ['similarity-analysis_from-FIR1_' cfg.sim_name]);
    end
end
% /sub-99 will be added below
% output of target folder will not be standard TDT content, because files can become too large

% add subject subdir
sourcedir = fullfile(sourcedir, sprintf('sub-%02i', cfg.subj)); % e.g. /derivatives/level1-FIR1/(func/)sub-99
targetdir = fullfile(targetdir, sprintf('sub-%02i', cfg.subj)); % e.g. /derivatives/similarityanalysis-timeresolved/sub-99

% add func subdir if it exists (for empirical data)
if isfolder(fullfile(sourcedir, 'func'))
    dispv(2, '  Found func directory, using it')
    sourcedir = fullfile(sourcedir, 'func'); % e.g. /derivatives/level1-FIR1/func
end

% check if folders exist (create target folder if not
assert(isfolder(sourcedir), 'Input directory not found: %s', sourcedir);
if ~isfolder(targetdir)
    [s, m] = mkdir(targetdir);
    assert(s, 'Could not create target folder %s, aborting. mkdir message: %s', targetdir, m);
end
disp(sourcedir); disp(targetdir); % show

%% Set input and output to cfg and load masks

fprintf('2 - Reading mask'); layout_line_break(1);
switch cfg.data_type
    case 'simulated'
        
        % Load full model mask to get whole voxel space dimensions
        mask_name = sprintf('sub-%02i_mask_fullmodel_%i.nii', cfg.subj, layer);
        maskdir = fullfile(sourcedir, mask_name);
        mask_hdr = spm_vol(maskdir);
        mask_vol = spm_read_vols(mask_hdr);
        
    case 'ROI'

        % Load ROI mask
        if isfield(cfg, 'ROI_dir')
            ROI_dir = cfg.ROI_dir;
        elseif isfield(cfg, 'corr_dir_ROI')
            [~, ROI_name, ~] = fileparts(cfg.corr_dir_ROI);
            if ~endsWith(ROI_name, '.nii'), ROI_name = [ROI_name, '.nii']; end
            ROI_dir = fullfile(cfg.preproc_dir, cfg.ROI_folder, ...
                sprintf('sub-%02i', cfg.subj), 'subjectspace_ROIs', ROI_name);
            if ~isfile(ROI_dir)
                warning('ROI mask %s missing! Skipping.', ROI_name); return; 
            end
        else
            warning('ROI mask missing! Skipping.'); return; 
        end
        mask_hdr   = spm_vol(ROI_dir);
        mask_vol   = spm_read_vols(mask_hdr);
        
    case 'wholebrain' 
        
        % Load standard wholebrain mask (usually called mask)
        mask_hdr = spm_vol(fullfile(sourcedir, 'mask.nii'));
        mask_vol = spm_read_vols(mask_hdr);
        
    otherwise
        error('unknown datatype for masks')
end
fprintf('    Reading mask done'); layout_line_break(2);

%% Set input and output to cfg: betas

% get betas in folder for verification of regressor information below
dispv(1, '1. Getting betas in folder...');
P = cellstr(spm_select('FPList', sourcedir, '^beta_[0-9]*.nii'));
assert(~isempty(P), 'No beta images found in %s, aborting', sourcedir);
hdrs = spm_vol(P);
fprintf('    Loading beta descriptions done.'); layout_line_break(2);
descrips = cellfun(@(x)x.descrip, hdrs, 'UniformOutput', false);

%% Load beta volumes

fprintf('3 - Reading %i volumes', length(hdrs)); layout_line_break(1);
vols = cell(length(hdrs), 1);
for h_ind = 1:length(hdrs)
    vols{h_ind} = spm_read_vols(hdrs{h_ind}, mask_vol);
end
fprintf('    Reading %i volumes done', length(hdrs));layout_line_break(2);

%% Keeping only voxels within mask
% Getting voxels within mask
% Returns a nimgs x 1 cell of nvxls x 1 data (nvxls: number of
% voxels in mask)

fprintf('4 - Masking volumes'); layout_line_break(1);
masked_vols = cellfun(@(x)x(mask_vol>0), vols, 'UniformOutput', false);

% Convert into one matrix nvxls x nimgs
masked_vols = cell2mat(masked_vols');
fprintf('    Masking volumes done.'); layout_line_break(2);

% Get rid of voxels with nans
masked_vols(any(isnan(masked_vols), 2), :) = [];

%% Delete trials in excess
% delete "Rest 36" trials. Ignores that the first sub has
% fewer images corresponding to certain trials. Also ignores that
% four participants have done only 7 runs.

% Get index of "~/Rest 36" trials
fprintf('5 - Getting indexes of trials to delete'); layout_line_break(1);
trials_to_filter = cellfun(@(x) regexp(x, 'Rest 36', 'match'), descrips, 'UniformOutput', false);
t_filter = cellfun(@isempty, trials_to_filter);
masked_vols_filtered = masked_vols(:, t_filter);
descrips = descrips(t_filter);
fprintf('    Data cleaned up.'); layout_line_break(2);

%% Throw an error if # voxels does not correspond to # voxels before SPM masking
% Make sure spm_global didn't exclude any voxel
% Store number of voxels

n_voxels = sum(mask_vol(:));

% fprintf('Expected/actual voxels: %i/%i\n', cfg.general_info.n_voxels_before_spm_masking, n_voxels);
% if n_voxels < cfg.general_info.n_voxels_before_spm_masking
%     fprinf('SUBJECT %i\n Mask name: %s', subj, maskname);
%     error('ATTENTION: SPM masked out some voxels. Check make_SPM_compatible function')
% elseif n_voxels > cfg.general_info.n_voxels_before_spm_masking
%     warning('ATTENTION: There are more voxels than expected! Check make_SPM_compatible function')
% end

%% 1 - Create list of descriptions divided by run
fprintf('6 - Create list of descriptions divided by run'); layout_line_break(1);

% Find out wether participants did 7 or 8 runs
check_eighth_run = cellfun(@(x) regexp(x, 'Sn\([8]\))', 'once','match'), descrips, 'UniformOutput', false);
if sum(~cellfun(@isempty, check_eighth_run))
    n_runs = 7;
else
    n_runs = 8;
end

list_runs = [];
descrips_list_runs = [];
for run = 1:n_runs
    % Output is 1 if the regex is included in the string, 0 if not
    list_runs{run} = ~cellfun(@isempty, regexp(descrips, ['Sn\(' num2str(run) '\) [^(constant)]'])); %#ok<*AGROW> 
    descrips_list_runs{run} = descrips(list_runs{run})';
    fprintf('    Run %i: %i scans.\n', run, size(list_runs{run},1)); layout_line_break(1);
end

%% Find indices of all even and odd runs

fprintf('7 - Find indices of all even and odd runs'); layout_line_break(1);
o_runs = [1 3 5 7];
if n_runs == 8
    e_runs = [2 4 6 8];
else
    e_runs = [2 4 6];
end

odd_runs = [];
for r= 1:length(o_runs)
    odd_runs{r} = list_runs{o_runs(r)};
end

even_runs = [];
for r= 1:length(e_runs)
    even_runs{r} = list_runs{e_runs(r)};
end

odd_idx = or(or(or(odd_runs{1}, odd_runs{2}), odd_runs{3}), odd_runs{4}); % extremely ugly, but it works
if n_runs == 8
    even_idx = or(or(or(even_runs{1}, even_runs{2}), even_runs{3}), even_runs{4});
else
    even_idx = or(or(even_runs{1}, even_runs{2}), even_runs{3});
end
odd_descr = descrips(odd_idx);
even_descr = descrips(even_idx);
% descrips_idx = or(odd_idx, even_idx);
fprintf('    Find indices of all even and odd runs done'); layout_line_break(2);

%% Use them here
fprintf('8 - Computing correlation matrix'); layout_line_break(1);

% correlation between all iamges, result nimgs x nimgs
switch  cfg.data_type % CHECK IF IT ACTUALLY MAKES A DIFFERENCE
    case 'ROIv'
        % dim of RHOs: nr of odd x even volumes   % TODO: CHANGE THIS
        RHOs = corr(masked_vols_filtered(:, odd_idx), masked_vols_filtered(:, even_idx ), 'rows', 'complete'); 
    otherwise
        % dim of RHOs: nr of odd x even volumes   % TODO: CHANGE THIS
        RHOs = corr(masked_vols_filtered(:, odd_idx), masked_vols_filtered(:, even_idx )); 
end

% Store number
fprintf('    Correlation matrix computed'); layout_line_break(2);

%% Extract condition names and calculating indeces using unique function
fprintf('9 - Extracting/cleaning up condition names and calculate position indices'); layout_line_break(1);

con_by_name_odd = cellfun(@(x) regexp(x, '[a-zA-Z]+[ ]*[0-9]{2,}*', 'match'), odd_descr);
con_by_name_even = cellfun(@(x) regexp(x, '[a-zA-Z]+[ ]*[0-9]{2,}*', 'match'), even_descr);

% Rename trial conditions so that they appear before their correspondent
% resting state
con_by_name_odd = strrep(con_by_name_odd, 'Trial', '.Trial');
con_by_name_even = strrep(con_by_name_even, 'Trial', '.Trial');
con_by_name_odd = strrep(con_by_name_odd, 'Cue', '.Cue');
con_by_name_even = strrep(con_by_name_even, 'Cue', '.Cue');

% Use "unique" function to: 1. Associate indices to trials; trials with the
% same name will have the same index 2. Get list of all trial names taken
% once
[con_unique_odd, ~, idx_unique_odd] = unique(con_by_name_odd(:));
[con_unique_even, ~, idx_unique_even] = unique(con_by_name_even(:));
fprintf('    Extracting/cleaning up condition names and calculate position indices done'); layout_line_break(2);

%% Compute average of correlation coefficients in RHOs correspondent to the same trials for simple/complex and odd or even runs;

% store the averages in a matrix
fprintf('10 - Average correlation coefficients in RHOs correspondent to the same trial within odd or even runs'); layout_line_break(1);
[cond_merged, odd_labelnames, even_labelnames] = average_mat(atanh(RHOs), ...
    idx_unique_odd, con_by_name_odd, idx_unique_even, con_by_name_even, @nanmean);
cond_merged = tanh(cond_merged);
fprintf('    Average matrix computed'); layout_line_break(2);

%% Store this together with the description

similarity_analysis_results.target_directory            = targetdir;
similarity_analysis_results.layer                       = layer;
similarity_analysis_results.n_voxels                    = n_voxels;
similarity_analysis_results.labelnames.odd              = odd_labelnames;
similarity_analysis_results.labelnames.even             = even_labelnames;
similarity_analysis_results.condition_unique.odd        = con_unique_odd;
similarity_analysis_results.condition_unique.even       = con_unique_even;
similarity_analysis_results.corr_results.RHOs           = RHOs;
similarity_analysis_results.corr_results.descrip_y      = odd_descr;
similarity_analysis_results.corr_results.descrip_x      = even_descr;
similarity_analysis_results.corr_results.datestr        = datestr(now);
similarity_analysis_results.corr_results.cond_merged    = cond_merged;

end
