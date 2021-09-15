function [X, L] = saa_get_des_mat(names, onsets, durations, pmod, settings, R)
% Create time courses for saa fMRI-like data
% FORMAT [X, L] = saa_get_des_mat(names, onsets, durations, pmod, settings, R)
%
%         NOTE: This function has been adapted from Joram Soch, taken from
%         https://github.com/JoramSoch/ITEM/blob/master/ITEM_get_des_mat.m
%         with permission from the author (ver 2021Mar12).
%         Modification: Added this and the above paragraph in the header.  
%         As orth was not implemented in Jorams version, an error was added.
%         The default for convolution is changed to trying to 
%         set settings.conv = 'stand'; because Joram suggested this. For 
%         convlution, only code for pmod has been kept, and it has data
%         demeaning has been taken out to keep original values (in pmod)
%         as hight of new data.
%
% Note that this function does NOT create parametric modulators, but uses
% the values in pmod (that are normally used to create parametric
% modulators) to create time series that use the pmod values as height of
% the to-be-convolved function % (or directly as value for setting.conv='none').
%
% INPUT
% 
%     names     - a  1 x c cell array of condition names
%     onsets    - a  1 x c cell array of condition onsets
%     durations - a  1 x c cell array of condition durations
%     pmod      - a  1 x c structure with the following fields:
%     o name    - a  1 x p cell array of modulator names
%     o param   - a  1 x p cell array of modulator values
%     note: orth has been removed because it is not used. would have been
%               - a  1 x c logical vector indicating orthogonalization
%                 note: currently not used, only '' is allowed as input
%     settings  - a structure variable with the following fields:
%     R         - an n x r design matrix of multiple regressors (e.g. RPs)
%                 Note: Changed order of settings and R
%     o n       - an integer, the number of scans
%     o TR      - a  scalar, the fMRI repetition time
%     o dt      - a  scalar, the microtime resolution
%     o hrf     - a  string, the convolution function (e.g. 'spm_hrf')
%     o conv    - a  string, the convolution mode ('none'/'stand')
%                 note: 'persist' from Jorams code removed
%     o mc      - a  logical, indicating mean-centering (for R)
%     o RPs     - a  logical, indicating realignment parameters (for R)
% 
%     X         - an n x p first-level fMRI design matrix
%     L         - a  1 x p cell array of regressor labels
% 
% FORMAT [X, L] = ITEM_get_des_mat(...) takes names, onsets, durations,
% pmod, orth and R in well-known SPM format and creates a design matrix
% following the specified settings (see below).
% 
% Adapted by Kai Goergen, with permission from Joram Soch
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% Modified: Kai, 2021/03/12
% First edit: 14/11/2017, 17:35 (V0.0)
%  Last edit: 22/11/2018, 11:45 (V0.1)


%% Set default values
%-------------------------------------------------------------------------%
% Note: orth removed. if your want to use it, put it e.g. at the end
% if isempty(orth) || nargin < 5              % orthogonalization
%     % orth = num2cell(false(size(names))); % original code from SPM, not yet used here    
% else
%     error('orth currently needs to be empty when provided, which means no orthogonalization will be performed. orthogonalization is not implemented in this version of the function yet. please check the original SPM function how to added it.')
% end
if isempty(R) || nargin < 5                 % additional regressors
    R = [];
end
if isempty(settings) || nargin < 6          % specification settings
    settings = struct([]);
end

%% Enact default settings
%-------------------------------------------------------------------------%
if ~isempty(settings)
    
    % Timing information
    %---------------------------------------------------------------------%
    % change from Kai:  number of images to create (n) and TR should be 
    % passed to function, removed default
%     if ~isfield(settings,'n'),    settings.n    = 1000;      end
%     if ~isfield(settings,'TR'),   settings.TR   = 2;         end
    
    % Onset regressors & parametric modulators
    %---------------------------------------------------------------------%
    if ~isfield(settings,'dt'),   settings.dt   = 0.01;      end
    if ~isfield(settings,'hrf'),  settings.hrf  = 'spm_hrf'; end
    if ~isfield(settings,'conv'), settings.conv = 'stand';   end
    
    % Additional regressors
    %---------------------------------------------------------------------%
    if ~isfield(settings,'mc'),   settings.mc   = true;      end
    if ~isfield(settings,'RPs')
        if size(R,2) == 6,        settings.RPs  = true;      end
        if size(R,2) ~= 6,        settings.RPs  = false;     end
    end
    
end

%% Specify time vector
%-------------------------------------------------------------------------%
t = 0:settings.dt:(settings.n*settings.TR);
z = zeros(1,2*numel(t)-1);
TRdt = round(settings.TR/settings.dt);

% Specify, normalize and partition HRF
%-------------------------------------------------------------------------%
if strcmp(settings.hrf,'spm_hrf')
    HRF = spm_hrf(settings.dt)';
else
    eval(strcat('HRF = ',settings.hrf,'(settings.dt);'));
end
 HRF     = HRF./max(HRF);

%% Preallocate design and labels
%-------------------------------------------------------------------------%
X = [];
L = [];

%% Create design matrix
%-------------------------------------------------------------------------%
for i = 1:numel(names)
    
    % Onset regressor
    %---------------------------------------------------------------------%
    o = round(onsets{i}./settings.dt);
    d = round(durations{i}./settings.dt);
    
    % the standard convolution is not necessary, we are only interested in
    % modulated regressors. 
    % NOTE: the stand / none part is different below, consider adapting it
    % if you uncomment these lines
%     % "standard" convolution
%     if strcmp(settings.conv,'stand') || strcmp(settings.conv,'none')
%         y = z; % init current conv vector y with empty vector z
%         for k = 1:numel(o)
%             y((o(k)+1):(o(k)+d(k))) = 1;
%         end
%         x = conv(y,HRF);
%         x = x(1:numel(z));
%     end
%     % no convolution
%     if strcmp(settings.conv,'none')
%         x = y;
%     end
%     % add regressor
%     X = [X, x(1:TRdt:(settings.n-1)*TRdt+1)'];
%     L = [L, names(i)];
    
    % Check for PMs
    %---------------------------------------------------------------------%
    if i <= numel(pmod)
        if ~isempty(pmod(i).name)
            for j = 1:numel(pmod(i).name)

                % Parametric modulator
                %---------------------------------------------------------%
                p = pmod(i).param{j};
%                 p = p-mean(p); % IMPORTANT MODIFICATION: dont use demeaning for data creation
                % create model time series (before convolution)
                y = z;
                for k = 1:numel(o)
                    y((o(k)+1):(o(k)+d(k))) = p(k);
                end
                % "standard" convolution
                if strcmp(settings.conv,'stand')                    
                    x = conv(y,HRF);
                    x = x(1:numel(z));
                % no convolution
                elseif strcmp(settings.conv,'none')
                    x = y;
                else
                    error('Unknown convolution methods settings.conv = %s', settings.conv)
                end
                % add modulator
                X = [X, x(1:TRdt:(settings.n-1)*TRdt+1)'];
                L = [L, {sprintf('%s x %s',names{i}, pmod(i).name{j})}];
            end
        end
    end
end

%% Add additional regressors
%-------------------------------------------------------------------------%
if ~isempty(R)
    if settings.mc
        R = R - repmat(mean(R),[settings.n 1]);
    end
    if settings.RPs
        R(:,end-2:end) = R(:,end-2:end) * (180/pi);
    end
    for i = 1:size(R,2)
        X = [X, R(:,i)];
        L = [L, {strcat('R',num2str(i))}];
    end
end