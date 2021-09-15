function [timeresolved_data, labelnames, ext_labelnames, masks, saa_cfg] = afx_saatable2timecourses(saa_cfg, datatable, TR, betas_per_sess)

% function [timeresolved_data,labelnames, ext_labelnames,  masks, saa_cfg] = afx_saatable2timecourses(saa_cfg, datatable, TR, betas_per_sess)
% IN
%   saa_cfg: struct with saa_cfg.testcase_definition
%   datatable: saa_datatable
%   betas_per_sess: n_sess x 1 vector, with number of betas per session. if
%       all sessions have the same number of images, one value is enough.
%   TR: scalar, fMRI repetition time in s
% OUT
%   timeresolved_data{sess_ind}: array with nsess matrices of size
%       betas_per_sess(sess_ind) x numel(labelnames)
%   labelnames: cellstr with one name per column of for matrix in 
%       timeresolved_data, as in saa_cfg.testcasedef, 
%       e.g. {'cond1', 'cond2', ...}
%   ext_labelnames{sess_ind}: cellstr with one name per column of for matrix in 
%       timeresolved_data{sess_ind}, added info which value was used to
%       define the height of each stick function during convolvution (like 
%       pmod in SPM), typically the same, e.g. {'cond1 x cond1', ...}
%   masks: struct array with one mask per testcase from saa_cfg.testcasedef
%       .maskname: testcasename
%       .mask_idx: indices of mask, e.g. [5 7] if the 5th and 7th column in
%           timeresolved_data contain the data for the current mask
%       .descrip: description of this testcase
%       .expect: expectation of this testcase
%       .labelnames: cellstr with used labelnames
%       .ext_labelnames: cellstr with used ext_labelnames
%   saa_cfg: saa_cfg with added field saa_cfg.map_saatable_values, with:
%       saa_cfg.map_saatable_values.sess(sess_ind).settings for 
%       saa_get_des_mat for each session
%
% Kai, 2021-03-12

%% get unique testentries
testcasedef = saa_cfg.testcase_definition;
unique_TableVariableNames = {};
for row_ind = 1:size(testcasedef, 1)
    curr_TableVariableNames = testcasedef.TableVariableNames{row_ind};
    % maybe todo: make sure there is no extra SAAVariableName left
    unique_TableVariableNames = unique([unique_TableVariableNames, curr_TableVariableNames]);
end

%% get timecourses for unique entries
if isfield(saa_cfg, 'session_var')
    sess = datatable.(saa_cfg.session_var);
else
    warning('cfg.session_var not defined. assuming all data comes from the same session')
    sess = ones(size(datatable, 1), 1); % no explicit session variable, put all in one
end
unique_sess = unique(datatable.(saa_cfg.session_var));

if numel(betas_per_sess) == 1 
    betas_per_sess = repmat(betas_per_sess, 1, numel(unique_sess));
end

for sess_ind = 1:length(unique_sess)
    curr_sess = unique_sess(sess_ind);
    curr_dat = datatable(sess == curr_sess, :);
    
    % special for this data: onset is stored slightly strange
    % onset is only onset for trial
    if ~isnumeric(curr_dat.onset)
        curr_dat.onset = str2double(curr_dat.onset);
    end
    if any(isnan(curr_dat.onset))
        error('TODO: REMOVE for general SAA: found nans, probably because rest data is used and this needs another onset field');
        error('Found nans in onsets, should not happen here, please check')
    end
    
    % init variables for this session
    clear names onsets durations pmod orth R
    for tc_ind = 1:numel(unique_TableVariableNames)
        curr_col_name = unique_TableVariableNames{tc_ind};
        curr_col_dat = curr_dat.(curr_col_name);
        % define info for regressor name, etc 
        % 'spmregressorname'can actually equal 'column'

        % question: how to use randn data in timecourses? what should each randn
        % modulate? 
        % I guess best would be height of the function
        % - i.e. parametric modulation in spm 
        % in a batch onset.mat, this can be stored in the struct array variable 

        if ~isnumeric(curr_col_dat)
            % try to get numeric data from string
            curr_col_dat_num = str2double(curr_col_dat);
            if ~any(isnan(curr_col_dat_num))
                % could convert numeric data, taking this
                curr_col_dat = curr_col_dat_num;
            end
        end
        if ~isnumeric(curr_col_dat)
            if ~exist('alreadywarned', 'var'), alreadywarned = {}; end
            if ~any(strcmp(curr_col_name, alreadywarned))
                % using unique to map string
                fprintf('afx_saatable2timecourses: Data in %s seem to be strings: using unique() to map to numbers\n', curr_col_name);
                % enough to warn once
                alreadywarned{end+1} = curr_col_name; %#ok<AGROW>
            end    
            [~, ~, curr_col_dat] = unique(curr_col_dat);
        end
            
        % specify information to create time series
        labelnames{tc_ind}    = unique_TableVariableNames{tc_ind};
        onsets{tc_ind}        = curr_dat.onset; % always the same at the moment, onset of trial
        durations{tc_ind}     = curr_dat.duration; % always the same at the moment, duration of trial (what to do with rest?)
        pmod(tc_ind).param{1} = curr_col_dat; % height of hrf, e.g. for binary variables 0 (no hrf) or 1 (standard height)
        pmod(tc_ind).name{1}  = unique_TableVariableNames{tc_ind}; % 
            % - pmod with pmod(onsetind).param{1}=RT pmod(onsetind).poly{1}=1 
            %   pmod(onsetind).name{1}='RT'
        R = []; % could carry e.g. realignment parameters
    end
    %% get all designs for this session
    settings = [];
    settings.n  = betas_per_sess(curr_sess);  % n       - an integer, the number of scans
    settings.TR = TR;   % TR      - a  scalar, the fMRI repetition time in s
    settings.dt = 1/33; % dt      - a  scalar, the microtime resolution
%     settings.conv = 'none';
    % if not set, defaults for:
    %     o hrf     - a  string, the convolution function (e.g. 'spm_hrf')
    %     o conv    - a  string, the convolution mode ('none'/'stand')
    % not used: stuff for R (that can e.g. take realignment parameters)
    %     o mc      - a  logical, indicating mean-centering (for R)
    %     o RPs     - a  logical, indicating realignment parameters (for R)

    [X, L] = saa_get_des_mat(labelnames, onsets, durations, pmod, settings, R);
    % normalise mean 0, std 1 (or all 0 if the values are all the same, i.e. const)
    % X = zscore(X);
    
%     % visualise (if you like)
%     %% plot (maybe with color limits)
%     figure('name', sprintf('sess %i (color limits)', curr_sess), 'position', [1           1        1600         800])
%     pX = X;
%     pL = L;
%     imagesc(pX);
%     title(sprintf('sess %i', curr_sess));
%     ylabel(sprintf('beta image (%i)', size(pX, 1)));
%     set(gca, 'XTick', 1:size(pX, 2), 'XTickLabel', pL, 'XTickLabelRotation', 90);
%     set(findall(gcf,'-property','FontSize'),'FontSize', 6)
%     set(gca, 'clim', quantile(pX(:), [.3 .97]))
%     colorbar
%     drawnow
% 
%     
    %% collect info to return
    saa_cfg.map_saatable_values.sess(sess_ind).sess = curr_sess;
    saa_cfg.map_saatable_values.sess(sess_ind).settings = settings;
    
    timeresolved_data{sess_ind} = X; %#ok<AGROW>
    
    if ~exist('ext_labelnames', 'var')
        ext_labelnames = L;
    else
        assert(isequal(ext_labelnames, L), 'Labelnames differ between session 1 & session %i (sess_ind %i), please check', curr_sess, sess_ind);
    end
end

%% Create masks
testcasedef = saa_cfg.testcase_definition;
masks = struct;
for row_ind = 1:size(testcasedef, 1)
    curr_testcasename = testcasedef.TestCaseName{row_ind};
    curr_TableVariableNames = testcasedef.TableVariableNames{row_ind};
    try
        curr_descrips = [testcasedef.Description{row_ind}{:}];
    catch
        curr_descrips = '';
    end
    try
        curr_expects = [testcasedef.Expectations{row_ind}{:}];
    catch
        curr_expects = '';
    end
        
    if isempty(curr_TableVariableNames)
        error('Testcase %s: .TableVariableNames is empty, but should be provided at this position. If only .SaaVariableNames exist, please resolve this at the top of this function to get TableVariableNames', curr_testcasename);
    end
    
    % find currently matching columns
    [c1,mask_idx,ib1] = intersect(labelnames, curr_TableVariableNames);
    % make sure all entries in curr_TableVariableNames exist in labelnames
    assert(length(ib1) == length(unique(curr_TableVariableNames)), 'Testcase %s: Some entries in curr_TableVariableNames are not included in labelnames. That should not happen, please check.', curr_testcasename)
    
    masks(row_ind).maskname = curr_testcasename;
    masks(row_ind).mask_idx = mask_idx;
    masks(row_ind).descrip  = curr_descrips;
    masks(row_ind).expect   = curr_expects;
    masks(row_ind).labelnames     = curr_TableVariableNames;
    masks(row_ind).ext_labelnames = ext_labelnames(mask_idx);
    masks(row_ind).info     = ['Mask created in afx_saatabke2timecourses' datestr(now)];
end