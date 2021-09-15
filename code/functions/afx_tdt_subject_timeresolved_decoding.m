% function [res_filenames, mask_data, cfg, results] = afx_tdt_subject_timeresolved_decoding(cfg)
%
% Handles full design creation, data filtering and subject-specific
% decoding for both simulated and empirical data. Specifically, it loads
% beta regressor caches, handles different decoding design options (e.g.
% pairwise or not, and a variety of labelling and cross-class decoding),
% reformats and plots accuracy and confusion matrix data, and saves all
% data.
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
%   res_filenames.resultfiles{mask_ind}: filenames of written resultsfiles
%       containing m_ind parts of mask_data{m_ind} as ACC_mask_data
%           .cfg: filenames of written cfg
%           .fig_fname: name of written figure
%   mask_data{m_ind}: struct including data for each mask (TODO: describe
%       fields)
%   cfg: including all files and design matrices
%   results: struct with all acccuracies/confusion matrices
%
% Kai and Ingmar, May 2020

function [res_fname, ACC_mask_data, cfg, results] = afx_tdt_subject_timeresolved_decoding(cfg)
%% Dev

do.plotting = 1;
do.reverse_file_order = 0; % Keep 0 unless testing for any confounds due to file order (very unlikely)

%% Source and target folder
dispv(1, 'Starting afx_tdt_subject_timeresolved_decoding for cfg.subj %i', cfg.subj)

sourcedir = cfg.level1_dir; % e.g. /analysis/share/corinna_kai/afx/derivatives/level1-FIR1
targetdir = cfg.level2_dir; % e.g. /analysis/share/corinna_kai/afx/derivatives/decoding-timeresolved
% /sub-99 will be added below
% output of target folder will not be standard TDT content, because files can become too large

% e.g. /analysis/share/corinna_kai/afx/derivatives/level1-FIR1/(func/)sub-99
sub_sourcedir = fullfile(sourcedir, sprintf('sub-%02i', cfg.subj));

% e.g. /analysis/share/corinna_kai/afx/derivatives/decoding-timeresolved/sub-99
sub_targetdir = fullfile(targetdir, sprintf('sub-%02i', cfg.subj));

% add func subdir if it exists (for empirical data)
if isfolder(fullfile(sub_sourcedir, 'func'))
    dispv(2, '  Found func directory, using it')
    % e.g. /analysis/share/corinna_kai/afx/derivatives/level1-FIR1/func
    sub_sourcedir = fullfile(sub_sourcedir, 'func');
end

% check if folders exist (create target folder if not
assert(isfolder(sub_sourcedir), 'Input directory not found: %s', sub_sourcedir);
if ~isfolder(sub_targetdir)
    [s, m] = mkdir(sub_targetdir);
    assert(s, 'Could not create target folder %s, aborting. mkdir message: %s', sub_targetdir, m);
end

dispv(1, 'Source dir.: %s', sub_sourcedir)
dispv(1, 'Output dir.: %s', sub_targetdir)

%% Save command window output in txt file

curr_date = strrep(datestr(now), ' ', '_');
try
    diaryfile = fullfile(sub_targetdir, sprintf('decoding_pipeline_output_%s.txt', curr_date));
    diary(diaryfile);
catch
    warning(['could not create diary ' diaryfile])
end

%% TDT Decoding settings

if exist('cfg', 'var')
    dispv(1, '  Updating passed cfg with TDT defaults...');
    cfg = decoding_defaults(cfg); % add TDT defaults
else
    cfg = decoding_defaults; %set defaults
end

if isfield(cfg, 'tdt_testmode'), cfg.testmode = cfg.tdt_testmode; end

cfg.results.dir = sub_targetdir; % results might not be stored with tdt, see below, but better be safe

%% Additional settings

% Option for potential speedup (not standard in TDT yet)
if strcmp(cfg.data_type, 'simulated')
    cfg.dont_read_headers_in_decoding_describe_data = 1;
end

% Allow unbalanced data
if strcmp(cfg.analysis_type, 'standard_cv') || strcmp(cfg.analysis_type, 'standard_difficulty_generalization')
    cfg.design.unbalanced_data = 'ok'; % multiclass
end

% Set input and output: result directory
cfg.results.write = 0;
% Disables TDT writing to own results folder; writing results currently NOT with TDT
% cfg.results.dir =  % would contain the output directory

%% Set input and output to cfg: masks
% goal: set cfg.files.masks
% read files from spm_mat, spm_reduced mat, or as fallback, from betas (or betas_cached)

% ---------------------------------------------------------- %
[cfg, layer_descr] = get_masks_to_decode(cfg, sub_sourcedir);
% ---------------------------------------------------------- %

%% Set input and output to cfg: betas

% ---------------------------------------------------------- %
[regressor_names, ~] = load_and_validate_betas_for_decoding(sub_sourcedir);
% ---------------------------------------------------------- %

%% Specify design labelling for training and testing

% For explanation of all possible combinations of labels and cross-class
% designs, including across-task generalizations, look in this function
% ----------------------------------------------------------------------%
all_cfgs = assign_design_labels(cfg, regressor_names, sub_sourcedir);
% ----------------------------------------------------------------------%

%% Check for confusion matrix: Change order in which data is loaded by TDT
% to avoid that "aeh" is loaded first because the betaxxx images have lower
% numbers

if do.reverse_file_order
    all_cfgs{1}.REVERSED_FILE_ORDER = 'Reversed file order to test if this has an effect on the confusion matrix (because maybe the order in which libsvm gets the data changes things';
    warning('Reversed file order to test if this has an effect on the confusion matrix (because maybe the order in which libsvm gets the data changes things')
    rev_fields = {'name', 'chunk', 'label', 'labelname', 'descr', 'set', 'xclass'};
    for rev_ind = 1:length(rev_fields)
        all_cfgs{1}.files.(rev_fields{rev_ind}) = all_cfgs{1}.files.(rev_fields{rev_ind})(end:-1:1);
    end
end

%% Loop across all_cfgs

total_results = {}; % init

for acfgs_ind = 1:length(all_cfgs)
    cfg = all_cfgs{acfgs_ind}; % Open current cfg in cell array
    
    %% Clean design
    
    if isfield(cfg.design, 'reduced_test_set') && cfg.design.reduced_test_set % TEST SET MODE
        dispv(1, 'A. Test mode active, removing all but 2 timepoints...')
        remove_filter = ~contains(cfg.files.descr, {'Trial00', 'Trial02'});
        warning('TEST MODE: ONLY TRIAL00 AND TRIAL02');
        cfg.design.make_set2conditions_map_str = 1; % descriptive strings per accuracy element
    else
        dispv(1, 'A. Removing "Rest 36" regressors from design')
        remove_filter = contains(cfg.files.descr, 'Rest 36');
    end
    
    fprintf('  Found %i regressors to remove. Removing these...\n', sum(remove_filter))
    disp(cfg.files.descr(remove_filter));
    cfg.files.name  = cfg.files.name(~remove_filter);
    cfg.files.chunk = cfg.files.chunk(~remove_filter);
    cfg.files.label = cfg.files.label(~remove_filter);
    cfg.files.descr = cfg.files.descr(~remove_filter);
    
    if ~any(strcmp({'standard_cv', 'standard_difficulty_generalization'}, cfg.analysis_type)) % xclass
        cfg.files.xclass = cfg.files.xclass(~remove_filter);
    end
    
    if isfield(cfg.files, 'labelname')
        cfg.files.labelname = cfg.files.labelname(~remove_filter);
    else
        warning('cfg.files.labelname does not exist, probably older version of TDT')
    end
    
    %% Change chunks from session to odd/even to cv between odd/even trials
    dispv(1, 'B. Change chunks to cv between odd/even blocks...')
    
    cfg.files.chunk = mod(cfg.files.chunk+1, 2)+1;
    
    %% Getting time information from description
    dispv(1, 'C. Getting time information from description')
    
    % In our data (example: 'Sn(8) KonCCue  00*bf(1)'), the timepoint is coded as 'Cue', 'Trial', or 'Rest' and the following two numbers
    % verify we found one timpoint per string
    cfg.files.timepointstr = regexp(cfg.files.descr, '(Cue|Trial|Rest)( )*[0-9][0-9]', 'match');
    n_found = cellfun(@length, cfg.files.timepointstr);
    assert(all(n_found==1), 'Found more or less timepoints in the descriptions of the following files: %s', ...
        num2str(find(n_found~=1)));
    
    % Take only cell as cellstr
    cfg.files.timepointstr = cellfun(@(x)x{1}, cfg.files.timepointstr, 'UniformOutput', false);
    
    % Assign number to each unique time
    % Get order that we want (first Cue, then Trial, then Rest)
    % Trick: replace Rest by Zest, then Cue, Trial, Zest in alphabetical order
    cfg.files.timepointstr  = strrep(cfg.files.timepointstr, 'Rest', 'Zest');
    [~, ~, uniqueid]        = unique(cfg.files.timepointstr); % add unique (now sorted) indices
    cfg.files.timepoint     = uniqueid;
    cfg.files.timepointstr  = strrep(cfg.files.timepointstr, 'Zest', 'Rest'); % Revert string replace
    
    %% Design: create time resolved design
    
    if ~any(strcmp({'standard_cv', 'standard_difficulty_generalization'}, cfg.analysis_type)) % xclass
        cfg.design.template_func = 'make_design_xclass_cv'; % for crossclass validation
    end
    if isfield(cfg.design, 'template_func') && strcmp(cfg.design.template_func, 'make_design_xclass_cv')
        dispv(1, 'D. Starting timeresolved design construction using "%s" labeling...', cfg.analysis_type); layout_line_break(1);
    else
        dispv(1, 'D. Starting timeresolved design construction for decoding...'); layout_line_break(1);
    end
    % ----------------------------------------------------------------------%
    cfg.design = make_design_timeresolved(cfg); % the design to create
    % ----------------------------------------------------------------------%
    
    %% Care about not too many figures from plot design
    % Plot only design for first subject

    if isfield(cfg, 'subj') && cfg.subs_todo(1) ~= cfg.subj
        disp('Plotting only design for sbj 1 to avoid to many figures')
        try % close & delete the template design
            close(cfg.design.fighandles.template_design);
            cfg.design = rmfield(cfg.design, 'fighandles');
            dispv(1, '    Closed and removed figure handles in cfg');
        catch
            dispv(1, '    Could not rmfield(cfg.design, "fighandles")'); close all; 
        end
        try cfg = rmfield(cfg, 'fighandles'); end %#ok<TRYNC> % try this too
        cfg.plot_design = 0; % don't plot in TDT
    end
    
    %% --------------------------- DECODING -------------------------------%%
    
    dispv(1, '  E. Calling decoding(cfg)...'); layout_line_break(1);
    % ----------------------------------------------------------------------%
    [results, cfg] = decoding(cfg); % passed data not necessary, equally fast and safer to ask TDT to do it
    % ----------------------------------------------------------------------%
    
    % Put back in cell array
    all_cfgs{acfgs_ind}         = cfg;
    total_results{acfgs_ind}    = results;  %#ok<*AGROW> % Combine results
    
    % Close handles and check directories
    if isfield(cfg.design, 'fighandles')
        try
            dispv(2, '  Closing plots...');
            close(cfg.design.fighandles.plot_design);
            close(cfg.design.fighandles.template_design);
        catch
            close all;
            dispv(2, '  Could not close figures, dont worry if none are open.')
        end
    end
    if isfield(cfg.results, 'dir') && ~strcmp(sub_targetdir, cfg.results.dir)
        warning('targetdir and cfg.results.dir are not the same anymore.')
        display(cfg.results.dir); display(sub_targetdir);
        disp('  Storing results to targetdir ...')
    end
end

%% Combine decoded cfgs and results

% QUESTION TO INGMAR: What does/should this do?
% Code below certainly does not do what it should, so removed it

% if length(all_cfgs) > 1
%     for c_ind = 2:length(all_cfgs) % start at 2nd
%         
%         % Leave only the design and files in remaining cfgs
%         curr_cfg          = all_cfgs{c_ind}; % init
%         curr_cfg.design   = all_cfgs{c_ind}.design;
%         curr_cfg.files    = all_cfgs{c_ind}.files;
%         curr_cfg.software = 'SPM12';
%         curr_cfg.cfg_combination = ['Warning: all redundant info was removed ', ...
%             'since first cfg == other cfgs except for design and files (possibly)'];
%         all_cfgs{c_ind} = curr_cfg; % overwrite cfg on that index
%     end
% end

%% Reshape ACC.ACCs so it is condition & timepoint train x test
dispv(1, '\n4. Reshaping accuracy results to condition & timepoint train x test...\n');

all_mask_ACCs = {};
for m_ind = 1:length(cfg.files.mask) % all masks
    ACC_mask_data = struct; % ACC_mask_data is later saved as .mat PER MASK, not as cells across masks
    % Init target images and get mask info
    ACC_mask_data.maskfile = cfg.files.mask{m_ind};
    [~, maskname, ~] = fileparts(ACC_mask_data.maskfile);
    
    % Determine (standard) pairwise decoding
    if ~isfield(cfg.design, 'make_pairwise'), cfg.design.make_pairwise = 0; end
    switch cfg.design.make_pairwise
        case true
            
            u_pairs = unique(cfg.design.set2conditions_map(:, 1));
            assert(length(u_pairs) == u_pairs(end), ... % check
                'Total pair count ~= vector length of (scalar) comparison labels');
            
            % Init counters
            set_u1 = length(unique(cfg.design.set2conditions_map(:, 2)));
            set_u2 = length(unique(cfg.design.set2conditions_map(:, 3)));
            set_c1 = 1; set_sum = (set_u1*set_u2); set_c2 = set_sum;
            
            % For standard pairwise decoding it's only 1 result set
            % with three labels per step when disregarding difficulty
            for p_ind = 1:length(u_pairs) % all class comparisons
                curr_results.accuracy.set = total_results{1}.accuracy.set(set_c1:set_c2); % get current pair
                % --------------------------------------------------------------- %
                ACCs{p_ind} = reshape_accuracies(all_cfgs{1}, curr_results, p_ind, m_ind, maskname); % See subfunc
                % --------------------------------------------------------------- %
                dispv(1, '  A. Constructed ACC matrix: comparison %0i/%0i "%s" for mask %0i/%0i: "%s"', ...
                    p_ind, length(u_pairs), ACCs{p_ind}.title, m_ind, length(all_cfgs{1}.files.mask), maskname);
                set_c1 = set_c1 + set_sum; set_c2 = set_c2 + set_sum; % update counters
            end % comparisons
            
        case false 
            
            % Standard_cv (multiclass): only has 1 (all labels versus each other)
            if any(strcmp({'standard_cv', 'standard_difficulty_generalization'}, cfg.analysis_type))
                u_comps = unique(cfg.design.set2conditions_map(:, 1));
                assert(length(u_comps) == u_comps(end), ... % check
                    'Total comparison count ~= vector length of (scalar) comparison labels');
                u_comps = length(u_comps);    
            % In the case of xclass total_results == 3 with 2 labels each per step
            elseif contains(cfg.analysis_type, {'xclass', 'cross'}) || ...
                    (isfield(cfg.design, 'make_pairwise') && cfg.design.make_pairwise == 1)
                u_comps = length(total_results);
            else
                error('timeresolved_decoding:unclearcomps', 'Cannot determine # of comparisons')
            end
            
            for u_ind = 1:u_comps % All class comparisons
                % --------------------------------------------------------------- %
                % in this case, the input argument comp_ind is not u_ind but 1 in all cases
                ACCs{u_ind} = reshape_accuracies(all_cfgs{u_ind}, total_results{u_ind}, 1, m_ind, maskname);
                % --------------------------------------------------------------- %
                dispv(1, '  A. Constructed ACC matrix: comparison %0i/%0i "%s" for mask %0i/%0i: "%s"', ...
                    u_ind, u_comps, ACCs{u_ind}.title, m_ind, length(all_cfgs{u_ind}.files.mask), maskname);
            end % comparisons
    end %switch
    
    % Plot decoding accuracy images
    if do.plotting
        fig_title = sprintf('Accuracies (ACCs) for sub-%02i for %s', cfg.subj, maskname);
        % ----------------------------------------------------------------------%
        acc_fighdl = plot_subject_accuracies(cfg, ACCs, fig_title, ACC_mask_data.maskfile, layer_descr);
        % ----------------------------------------------------------------------%
        dispv(1, '  B. Result plotted for mask %0i/%0i: "%s"', m_ind, length(cfg.files.mask), maskname);
    end
    all_mask_ACCs{m_ind} = ACCs; % store
    
    % Save plot
    if isfield(results, 'accuracy') && do.plotting
        acc_plot_fname = fullfile(sub_targetdir, sprintf('decoding-accuracies_%s', maskname)); % filename
        save_plot(acc_plot_fname, acc_fighdl, 1)
        close(acc_fighdl);
    end
end % masks
clear ACC_mask_data % avoid reuse (contains only the last mask_data)

%% Confusion matrix setup

if (isfield(cfg.design, 'make_pairwise') && cfg.design.make_pairwise == 1 && ...
        isfield(results, 'confusion_matrix')) || ...
        ~any(strcmpi(cfg.analysis_type, {'standard_cv', 'standard_difficulty_generalization'}))
    warning([ 'Pairwise CONFMAT figures discontinued. If needed, run multiclass ', ...
        'with cfg.results.output = {"confusion_matrix"}']);
    
    % To avoid problems later on
    try results = rmfield(results, 'confusion_matrix'); end
    cfg.results.output = {'accuracy'};
end

% Create CONFMATs ------------------------------------------ %
[all_CONFMATs, ~] = afx_subject_timeresolved_confmats(cfg, results);
% ---------------------------------------------------------- %

%% Save data
dispv(1, '5. Saving data...');

% Save, assumes all cfgs are comparable, since saving all separate cfgs will
% require more changes later on
cfg_fname = fullfile(sub_targetdir, 'cfg');
disp(['  Saving LAST cfg: ', cfg_fname]);
save(cfg_fname, 'cfg');

% Saving all cfgs
all_cfg_fname = fullfile(sub_targetdir, 'all_cfg');
disp(['  Saving ALL cfgs: ', all_cfg_fname]);
save(all_cfg_fname, 'all_cfgs');

% Get accuracy set (e.g. 1 for multiclass, 3 for xclass) for all masks
for m_ind = 1:length(all_mask_ACCs) % loop over masks
    ACC_mask_data = []; % avoid leak from previous looper or earlier in the script
    for a_ind = 1:length(all_mask_ACCs{m_ind}) % loop over accuracy sets (e.g. 1 for multiclass, 3 for xclass)
        
        [~, mfname] = fileparts(cfg.files.mask{m_ind});
        
        ACC_mask_data.maskfile = cfg.files.mask{m_ind};
        ACC_mask_data.info = ['generated by mfilename: ' mfilename '. ' ... % not the current mask name, but the current .m-file
            'cfg.time_log.creation_date: ' cfg.time_log.creation_date ...
            '.files, etc copied from tdt cfg.'];
        ACC_mask_data.m_ind = m_ind;
        
        % no sub-xx_ duplicates (if they were present in mask name)
        if startsWith(mfname, 'sub-')
            mfname = mfname(8:end);
        end
        
        % Collect info for current mask (name already above)
        ACC_mask_data.ACCs{a_ind} = all_mask_ACCs{m_ind}{a_ind};
    end
    
    % Save
    dispv(1, '  Saving ACCs for mask "%s"', ACC_mask_data.ACCs{1}.maskname);
    res_fname = fullfile(sub_targetdir, sprintf('res-ACC_sub-%02i_%s', cfg.subj, mfname));
    save(res_fname, 'ACC_mask_data');
end % masks
clear ACC_mask_data % avoid reuse (contains only the last mask_data)

% Store confusion matrix if present
if exist('all_CONFMATs', 'var') && ~isempty(all_CONFMATs)
    for m_ind = 1:length(cfg.files.mask) % no comparison loop, as it is not supported
        
        % No sub-xx_ duplicates (if they were present in mask name)
        [~, mfname] = fileparts(cfg.files.mask{m_ind}); % get mask name again
        if startsWith(mfname, 'sub-'), mfname = mfname(8:end); end
        
        % Save
        CONFMAT = all_CONFMATs{m_ind};
        confmat_fname = fullfile(sub_targetdir, sprintf('res-CONFMAT_sub-%02i_%s', cfg.subj, mfname));
        dispv(1, '  Saving corresponding CONFMAT:\n    %s...', confmat_fname);
        save(confmat_fname, 'CONFMAT')
    end % masks
end

% Done
disp('ACC decoding done.')

%% Save diary

diary off

end

