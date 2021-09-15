% reduce_SPM_mats_from_dirs(cfg)
%
% Function to reduce all SPM mats from all participants' directories
% Using TDT function "reduce_SPM_filesize.m" written by Kai.
%
% Leonardo, 2020

% Updated and added to pipeline by Ingmar, 21-07-21

function reduce_SPM_mats_from_dirs(cfg)

% For all subjects
for s_ind = 1:length(cfg.subs_todo)
    sub = cfg.subs_todo(s_ind);
    
    % Folder: eg .../derivatives-sim/level1-FIR1_sim-suppression/sub-01/
    spmdir_todelete = fullfile(cfg.derivative_dir, [cfg.level1 '_' cfg.sim_name], ...
        sprintf('sub-%02i', sub));
    
    if isfile(fullfile(spmdir_todelete, 'SPM.mat'))
        
        % Custom function (but does not delete old SPM.mat)
        reduce_SPM_filesize(spmdir_todelete, 1);
        
        % Delete old one
        % delete(spmdir_todelete)
        fprintf('\n    SPM.mat(s) reduced in filesize for sub-%02i!', sub)
    else
        fprintf('\n    SPM.mat not present or already reduced for sub-%02i!', sub)
    end
end