% function ts_data = make_time_series_model_from_cue(cfg, SPM, conditions, model, modelname, hrf_conv, noise_before_conv, noise_after_conv)
%
% Create time series from content of an existing SPM.mat.
%
% INPUT
% ----------------------------------------------------------------------
% conditions : single string or array with the names of the regressors
%              (e.g. conditions = {'left', 'right'})
% SPM : SPM mat as struct (from load SPM.mat)
% hrf_conv : 'hrf_conv' or 'no_hrf_conv'
% noise_before_conv & noise_after_conv: noise parameters set prior and
%                                       passed to this function
%
% OUTPUT
% ----------------------------------------------------------------------
% dat : [1 x runs] struct storing design matrices for each run
%     dat(s).X                = X;
%     dat(s).regr_names       = [U(idx).name];
%     dat(s).convolution_type = hrf_conv;
%     dat(s).curr_run         = s;

% Function written using 'spm_fMRI_design' as reference -- by Leonardo
% Adapted Kai, 2020-02-25

function ts_data = make_time_series_model_from_cue(cfg, SPM, conditions, model, modelname, hrf_conv, noise_before_conv, noise_after_conv)
%% Adjust inputs

% Check if model is column vector
if size(model,1) == 1
    model = model';
end

% Get condition name from SPM mat
conditions_short_names = conditions; % store input conditions 

% Get current hypothesis name
sim_name = cfg.sim_name;

%% -Microtime onset and microtime resolution

try
    fMRI_T     = SPM.xBF.T;
    fMRI_T0    = SPM.xBF.T0;
catch
    fMRI_T     = spm_get_defaults('stats.fmri.t');
    fMRI_T0    = spm_get_defaults('stats.fmri.t0');
    SPM.xBF.T  = fMRI_T;
    SPM.xBF.T0 = fMRI_T0;
end

%% Time units, dt = time bin {secs} mode

SPM.xBF.dt  = SPM.xY.RT/SPM.xBF.T;
SPM.xBF     = spm_get_bf(SPM.xBF);

%% Get session specific design parameters

Xx      = [];
Xname   = {};
ts_data = struct();

%% Loop over all required time series

for ts_ind = 1:size(SPM.nscan, 2)
    
    % If model is the modelrest for the sticknoaftereffect_randstickrest
    % hypothesis, compute a new model (we want the sticks in the resting state to be randomised across participants and
    % runs)
    if contains(modelname , 'rest') && contains(cfg.sim_name, 'randstickrest')
        [~, ~, model] = models_from_hypothesis('sticknoaftereffect_randstickrest', cfg.saa.data_gen.dt, cfg.saa.data_gen.response, 0); % get only modelrest
    end 
    
    %-Number of scans for this session
    k = SPM.nscan(ts_ind);

    % Group onsets and causes by condition (cue + trial + resting state)
    % Update SPM.mat
    U = spm_get_ons(SPM, ts_ind);
    
    % Get name of the target condition
    names = [U.name];
    
    % Get full name of conditions 
    conditions = condition_to_model_from_cue(SPM, ts_ind, conditions_short_names); 

    % Find index conditions
    idx = ismember(names, conditions); 
    
    % When model is a cellstr
    if iscell(model); model = model{:}; end 
    
    % Add noise to the model if selected
    if exist('noise_before_conv', 'var') && all(noise_before_conv ~= 0)
        model = model + noise_before_conv*randn(size(model));
    end
    
    %% Convolve stimulus functions with basis functions
    
    switch hrf_conv
        case 'hrf_conv' % Convolve with hrf, if selected
            dispv(2, '    Convolving %0i/%0i time series with model and hrf', ts_ind, length(SPM.nscan))
            [X,Xn,Fc] = spm_Volterra(U(idx), SPM.xBF.bf, SPM.xBF.Volterra);
            if ~isempty(X)
                model_mat = [];
                for m_ind = 1:size(X,2)
                    model_mat(:, m_ind) = conv(model, X(:,m_ind)); %#ok<AGROW>
                end
                X = model_mat;
            end                      
        case 'no_hrf_conv' % No convolution
            dispv(2, '    Convolving %0i/%0i time series with model (no hrf)', ts_ind, length(SPM.nscan))
            [X,Xn,Fc] = spm_Volterra(U(idx),  model, SPM.xBF.Volterra);        
        otherwise
            error('Please indicate whether you want to convolve your model with the hrf or not.')
    end
    
    % Add random noise if selected 
    if exist('noise_after_conv', 'var') && noise_after_conv ~= 0
        X = X + noise_after_conv* randn(size(X));
    end

    %-Resample regressors at acquisition times (32 bin offset)
    if ~isempty(X)
        X = X((0:(k - 1)) * fMRI_T + fMRI_T0 + 32,:);
    end
     
    %-Get user specified regressors
    C     = SPM.Sess(ts_ind).C.C;
    Cname = SPM.Sess(ts_ind).C.name;
    
    %-Append mean-corrected regressors and names
    X     = [X spm_detrend(C)];
    Xn    = [Xn(:)' Cname(:)'];
    
    %-Session structure array
    SPM.Sess(ts_ind).U      = U;
    SPM.Sess(ts_ind).C.C    = C;
    SPM.Sess(ts_ind).C.name = Cname;
    SPM.Sess(ts_ind).col    = size(Xx,2) + (1:size(X,2));
    SPM.Sess(ts_ind).Fc     = Fc;
    
    %-Append into Xx and Xb
    Xx = blkdiag(Xx,X);
    
    %-Append names
    for i = 1:length(Xn)
        Xname{end + 1} = [sprintf('Sn(%i) ',ts_ind) Xn{i}];
    end

    % - Add to dat
    ts_data(ts_ind).X                = X;
    ts_data(ts_ind).regr_cue_names   = [U(idx).name];
    ts_data(ts_ind).convolution_type = hrf_conv;
    ts_data(ts_ind).curr_run         = ts_ind;
end

end

% Function to get onset of condition (that is the cue associated with it)
% as it is named in SPM.Sess.U.name

function cond_U_name = condition_to_model_from_cue(SPM, run, cond)

cond_U_name = cell(size(cond));

for c_ind = 1:length(cond)
    
    % Left or right
    if strcmp(cond{c_ind}, 'left') || strcmp(cond{c_ind}, 'right')
        cond_U_name{c_ind}  = cond{c_ind};
    else
        names = [SPM.Sess(run).U.name];
        curr_name = sprintf('cue-%s', cond{c_ind});
        cond_U_name{c_ind} = names(filter_condition(names, curr_name));
        cond_U_name{c_ind} = cond_U_name{c_ind}{:};
    end
end

end