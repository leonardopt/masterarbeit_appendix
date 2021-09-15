% all_paths = set_all_paths()
%
% Function to set all paths. Adds SPM and TDT to the path if available.
%
% Note: global scope for output might be removed in a future MATLAB
% release.
%
% Kai, 2019

function all_paths = set_all_paths()
%% Paths

global all_paths
if isempty(all_paths)
    global_all_paths();
end

%% Add TDT and SPM and saa-tdt

addpath(genpath(all_paths.tdt));
addpath(all_paths.spm12);
sprintf('SPM path: %s\n', all_paths.spm12);
try
    addpath(all_paths.tdtsaa);
    display(all_paths.tdtsaa)
catch
    warning('Set all_paths.tdtsaa in global_all_paths_private (and call global_all_paths_private manually once to update all_paths)')
end
%% Add folder from this repo

addpath(all_paths.svn);
addpath(genpath(all_paths.svn));

%% Add remote setup

if ispc
    addpath('z:\share\software\matlab_libraries\remote_setup')
else
    addpath('/analysis/share/software/matlab_libraries/remote_setup')
end

end