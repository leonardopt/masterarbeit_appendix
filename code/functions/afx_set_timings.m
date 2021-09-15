% [customticks, customlabels] = afx_set_timings(few, timepoint_total)
%
% Sets custom ticks and timings for decoding results plotting.
% 
% Note: cue is actually 3/4(?) seconds, but first trial TR locks to
% the first trial onset, which is 7.5s in the full recording. Since
% we only accounted for 1 cue regressor (Cue 00), this regressor
% actually only contains 1.5s of cue and 0.5s of the first
% pre-trial fixation cross.
%
% Ingmar, 09-2020

function [customticks, customlabels] = afx_set_timings(few, timepoint_total)

if ~exist('timepoint_total', 'var'), timepoint_total = 43; end

% Tick info
cue_len         = 2; % s
fixate          = 0.5; % s
trials_start    = cue_len + fixate;
trials_end      = cue_len + (6 + fixate) * 8; % s
rest_end        = timepoint_total * 2; % -1 to keep a tick at start of TR
trialtimes      = (trials_start:(6 + fixate):trials_end)/2; % t1-t2 start = 6.5s, divided by TR
resttimes       = round((trials_end:8:rest_end)/2);

% All or less ticks to display for more font space  
if ~exist('few', 'var') || isempty(few) || few == 0
    customticks     = [trialtimes resttimes];
    customlabels    = [compose('T%01i', 1:8), compose('R%is', ...
                        resttimes*2 - trials_end)];     
elseif few == 1   
    trialtimes      = trialtimes(2:2:6);
    resttimes       = linspace(trials_end, trials_end + 32, 3)/2;
    customticks     = [trialtimes resttimes];
    customlabels    = [compose('T%i', 2:2:6), compose('R%is', ...
                        resttimes*2 - trials_end)];
end