% sets all_paths as global variable depending on user
%
% Author: kai, version: 20-02-13
%
% IMPORTANT
%
% This file should only exist one time per project, in it's base folder. It
% should NOT be copied into subfolders.
%
% What should can be copied into each subfolder is the file
% global_all_paths.m. This then searchs for this _private.m file and adds
% the specified paths.
%
% This file contains a list of directories where information related to the
% project is found. This can be different for different users and different
% computers.
%
% If you are a new user or on a new computer, please check your username
% and hostname (machine name), and copy and modify the path to your need.
% 
% SEE ALSO: global_all_paths.m
%
% Kai, 20-02-13


function global_all_paths_private()

global all_paths

% get username & machine you are working on (in case you woprk on multiple machines)
if isunix
    % on mac, isunix == 1, so it works there, too
    [dummy, user_name] = system('whoami'); % exists on every unix that I know of
    [dummy, host_name] = system('hostname'); % exists on every unix that I know of
elseif ispc % pc means windows only
    [dummy, user_name] = system('echo %USERNAME%');
    [dummy, host_name] = system('echo %USERDOMAIN%');
else
    error('Unkown OS')
end

% remove empty characters at front & end
user_name = strtrim(user_name);
host_name = strtrim(host_name);

if strcmp(user_name, 'Kai-Scioi')
    all_paths.svn = 'C:\kai\afx_ingmar_leo\afx_timeresolved';
    all_paths.needed_functions = 'C:\kai\afx_ingmar_leo\TimeResolvedAFXSimSAA\needed_functions';
    all_paths.svn_cluster = '/home/kai/projects/afx_timeresolved';
    all_paths.localdir = 'C:\kai\';
    all_paths.tdt = 'C:\kai\tdt\decodingtoolbox-code\decoding_toolbox';
    all_paths.spm12 = 'C:\Users\Kai-Scioi\Documents\MATLAB\spm12';
    all_paths.data_base_dir = 'z:/share/corinna_kai/afx';
    all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx'; % where cluster look for the data
elseif strcmp(user_name, 'kai')
    if strcmp(host_name, 'Kai-Gorgens-MacBook.local')
        all_paths.svn = '/Users/kai/Documents/!Projekte/afx_timeresolved';
        all_paths.svn_cluster = '/home/kai/projects/afx_timeresolved';
        all_paths.localdir = '/local/kai';
        all_paths.tdt = '/Users/kai/Documents/!Projekte/Decoding_Toolbox/trunk/decoding_betaversion';
        all_paths.spm12 = '/Users/kai/Documents/!Studium/Software/matlab libraries/spm12';
        all_paths.data_base_dir = '$NOT_AVAILABLE: Set in global_all_paths_private$';
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx'; % where cluster look for the data
    elseif isunix
        all_paths.svn = '/home/kai/projects/afx_timeresolved'; % '/home/kai/projects/aeffect_218';
        all_paths.svn_cluster = '/home/kai/projects/afx_timeresolved'; % ONLY for condor job generation: path were the same file are on the cluster
                                                            % SETUP A
                                                            % DIFFERENT
                                                            % PROFILE FOR
                                                            % cluster, IF
                                                            % YOU USE IT!

        all_paths.localdir = '/scratch/condor_tmp/kai/afx';
        all_paths.tdt = '/home/kai/projects/decoding_tool/trunk/decoding_betaversion';
        all_paths.spm12 = '/analysis/share/spm12';
        all_paths.tdtsaa = '/home/kai/SAAdev/tdt-saa/TDT_SAA_dev/TDT_SAA';
        all_paths.data_base_dir = '/analysis/share/corinna_kai/afx';
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx/'; % where cluster look for the data
    else
        all_paths.svn = 'H:\projects\afx_timeresolved';
        all_paths.svn_cluster = '/home/kai/projects/afx_timeresolved'; % ONLY for condor job generation: path were the same file are on the cluster
                                                            % SETUP A
                                                            % DIFFERENT
                                                            % PROFILE FOR
                                                            % cluster, IF
                                                            % YOU USE IT!

        all_paths.localdir = 'C:\temp';
        all_paths.tdt = 'H:\projects\decoding_tool\trunk\decoding_toolbox\';
        all_paths.spm12 = 'z:\share\spm12\';
        all_paths.data_base_dir = 'z:\share\corinna_kai\afx\';
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx'; % where cluster look for the data
    end

    
    
elseif strcmp(user_name, 'anthony') || strcmp(user_name, 'anthony819')
    if strcmp(host_name, 'Anthonys-MacBook-Pro.local')
        all_paths.svn = '/Users/anthony819/Documents/Germany/BCCN Lab/aeffect_218/root';
        all_paths.svn_cluster = '/home/anthony/projects/aeffect_218';
        all_paths.localdir = '/Users/anthony819/Documents/MATLAB/BCCN Lab/temp';
        all_paths.tdt = '/Users/anthony819/Documents/MATLAB/decoding_toolbox_v3.994';
        all_paths.spm12 = '/Users/anthony819/Documents/MATLAB/spm12';
        all_paths.data_base_dir = '/Users/anthony819/Documents/MATLAB/BCCN Lab/aeffect';
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx'; % where cluster look for the data
    elseif isunix
        all_paths.svn = '/home/anthony/projects/aeffect_218';
        all_paths.svn_cluster = '/home/anthony/projects/aeffect_218'; % ONLY for condor job generation: path were the same file are on the cluster
                                                            % SETUP A
                                                            % DIFFERENT
                                                            % PROFILE FOR
                                                            % cluster, IF
                                                            % YOU USE IT!

        all_paths.localdir = '/scratch/condor_tmp/anthony/afx';
        %all_paths.localdir = '/home/anthony/temp/afx';
        all_paths.tdt = '/home/anthony/projects/decoding_toolbox_v3.996';
        all_paths.spm12 = '/analysis/share/spm12';
        all_paths.data_base_dir = '/analysis/share/corinna_kai/afx';
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx/'; % where cluster look for the data
    else
        all_paths.svn = 'H:\projects\aeffect_218';
        all_paths.svn_cluster = '/home/anthony/projects/aeffect_218'; % ONLY for condor job generation: path were the same file are on the cluster
                                                            % SETUP A
                                                            % DIFFERENT
                                                            % PROFILE FOR
                                                            % cluster, IF
                                                            % YOU USE IT!

        all_paths.localdir = 'C:\temp';
        all_paths.tdt = 'H:\projects\decoding_toolbox_v3.996\';
        all_paths.spm12 = 'z:\share\spm12\';
        all_paths.data_base_dir = 'z:\share\corinna_kai\afx\';
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx'; % where cluster look for the data
        
    end


elseif strcmp(user_name, 'corinna') || strcmp(user_name, 'corinnapehrs')
    if strcmp(host_name, 'Corinnas-MacBook-Pro.local')
        all_paths.svn = '/Users/corinnapehrs/Documents/projects/aeffects';
        all_paths.svn_cluster = '/home/corinna/aeffects';
        %all_paths.localdir = '/Users/anthony819/Documents/MATLAB/BCCN Lab/temp';
        all_paths.tdt = '/Users/corinnapehrs/Documents/projects/affects/decoding_toolbox_v3.997';
        all_paths.spm12 = '/Users/corinnapehrs/spm12';
        %all_paths.data_base_dir = '/Users/anthony819/Documents/MATLAB/BCCN Lab/aeffect';
        %all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx'; % where cluster look for the data
    elseif isunix
        all_paths.svn = '/home/corinna/aeffects';
        all_paths.svn_cluster = '/home/corinna/aeffects'; % ONLY for condor job generation: path were the same file are on the cluster
                                                            % SETUP A
                                                            % DIFFERENT
                                                            % PROFILE FOR
                                                            % cluster, IF
                                                            % YOU USE IT!

        all_paths.localdir = '/scratch/condor_tmp/corinna/afx';
        %all_paths.localdir = '/home/anthony/temp/afx';
        all_paths.tdt = '/home/corinna/aeffects/decoding_toolbox_v3.997';
        all_paths.spm12 = '/analysis/share/spm12';
        all_paths.data_base_dir = '/analysis/share/corinna_kai/afx';
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx/'; % where cluster look for the data
    else
        all_paths.svn = 'H:\projects\aeffect_218';
        all_paths.svn_cluster = '/home/anthony/projects/aeffect_218'; % ONLY for condor job generation: path were the same file are on the cluster
                                                            % SETUP A
                                                            % DIFFERENT
                                                            % PROFILE FOR
                                                            % cluster, IF
                                                            % YOU USE IT!

        all_paths.localdir = 'C:\temp';
        all_paths.tdt = 'H:\projects\decoding_toolbox_v3.996\';
        all_paths.spm12 = 'z:\share\spm12\';
        all_paths.data_base_dir = 'z:\share\corinna_kai\afx\';
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx'; % where cluster look for the data
    end

    
elseif strcmp(user_name, 'leonardo.pettini') || strcmp(user_name, 'pettinil')
    if strcmp(host_name, 'Leonardos-MBP.lan') || strcmp(host_name, 'Leonardos-MacBook-Pro.local')
        all_paths.svn = '/Users/leonardo.pettini/Documents/MATLAB/aeffect_218';
        % all_paths.svn_cluster = '/home/pettinil/';
%         all_paths.localdir = '/Users/anthony819/Documents/MATLAB/BCCN Lab/temp';
        all_paths.tdt = '/Users/leonardo.pettini/Documents/MATLAB/decoding_toolbox_v3.997';
        all_paths.spm12 = '/Users/leonardo.pettini/Documents/MATLAB/spm12';
        all_paths.data_base_dir = '/Volumes/analysis/share/corinna_kai/afx';
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx'; % where cluster look for the data
    elseif isunix && strcmp(user_name, 'pettinil')
        all_paths.svn = '/home/pettinil/afx';
        all_paths.svn_cluster = '/home/pettinil/afx'; % ONLY for condor job generation: path were the same file are on the cluster
                                                            % SETUP A
                                                            % DIFFERENT
                                                            % PROFILE FOR
                                                            % cluster, IF
                                                            % YOU USE IT!

        all_paths.localdir = '/scratch/condor_tmp/pettinil/afx';
       % all_paths.tdt = '/home/pettinil/Documents/MATLAB/decoding_toolbox_v3.997';
        %all_paths.spm12 = '/home/pettinil/afx/Analysis/Research_project_2_BCCN/Toolboxes/spm12/';
        all_paths.spm12 = '/analysis/share/spm12/';
        all_paths.data_base_dir = '/analysis/share/corinna_kai/afx';
        all_paths.tdt = '/home/pettinil/afx/Analysis/Research_project_2_BCCN/Toolboxes/decoding_toolbox_v3.999'
        all_paths.needed_functions = '/home/pettinil/afx/Analysis/TimeResolvedAFXSimSAA/needed_functions'
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx/'; % where cluster look for the data
    else
        all_paths.svn = 'H:\projects\aeffect_218';
        %all_paths.svn_cluster = '/home/anthony/projects/aeffect_218'; % ONLY for condor job generation: path were the same file are on the cluster
                                                            % SETUP A
                                                            % DIFFERENT
                                                            % PROFILE FOR
                                                            % cluster, IF
                                                            % YOU USE IT!

        all_paths.localdir = 'C:\temp';
        all_paths.tdt = 'H:\projects\decoding_toolbox_v3.996\';
        all_paths.spm12 = 'z:\share\spm12\';
        all_paths.data_base_dir = 'z:\share\corinna_kai\afx\';
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx'; % where cluster look for the data
        
    end
    

elseif strcmp(user_name, 'ingmar') || strcmp(user_name, 'bccnwin\ingmar')
    if strcmp(host_name, 'BCCNWIN')
        all_paths.svn = 'C:\LocalDir_Ingmar\Research_project_2_BCCN\';
        all_paths.tdt = 'C:\LocalDir_Ingmar\Toolboxes\decoding_toolbox_v3.999\';
        all_paths.tdt_saa = 'C:\LocalDir_Ingmar\Toolboxes\TDT_SAA_rev443_v0_9_Ingmar\';
        all_paths.spm12 = 'C:\LocalDir_Ingmar\Toolboxes\spm12\';
        all_paths.data_base_dir = 'Z:\share\corinna_kai\afx\';
        all_paths.data_base_dir_nathan = 'Z:\share\corinna_kai\afx\'; % 
        all_paths.localdir = 'C:\LocalDir_Ingmar\';
    elseif isunix
        all_paths.svn = '/home/ingmar/TimeResolvedAFXSimSAA';
        all_paths.svn_cluster = all_paths.svn; %only working on unix in cluster
        all_paths.tdt = '/home/ingmar/Toolboxes/tdt_3_999B_beta';
        all_paths.tdt_saa = '/home/ingmar/Toolboxes/TDT_SAA_rev443_v0_9_Ingmar';
        all_paths.spm12 = '/home/ingmar/Toolboxes/spm12';
        all_paths.data_base_dir = '/analysis/share/corinna_kai/afx';
        all_paths.data_base_dir_nathan = '/analysis/share/corinna_kai/afx';
        all_paths.localdir = '/home/ingmar/local_H_folder'; % "fake" local folder
    end
    
        
    
    
    
    
% elseif strcmp(user_name, 'NEWUSER')  

% copy allpaths.* from above and add
% your path here
    
    
else
    error('Username %s is unkown, please check and add to global_all_paths_private.m', user_name)
end

%% Check for unspecified defaults

% if ~isfield(all_paths, 'tmp_dir')
%     if ispc
%         all_paths.tmp_dir = fullfile('C:\tmp\', user_name);
%     else
%         all_paths.tmp_dir = fullfile('/tmp/', user_name);
%     end
%     warning('global_all_paths_private:all_paths_tmp_dir_unset', 'all_paths.tmp_dir has not been set explicitely in global_all_paths_private. Using default all_paths.tmp_dir = %s', all_paths.tmp_dir);
% end
        
%% Report all path used
disp('Global variable all_paths set as follows:')
display(all_paths)