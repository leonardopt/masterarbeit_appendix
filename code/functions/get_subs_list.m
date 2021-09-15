% subs_todo = get_subs_list()  
%
% Function to get list of subjects. Relies on "all_paths" and SPM.
% Optional input: indices of a subset of subjects.
%
% Kai, 2019

% Updated to pass subjects for exclusion - Ingmar, March 21

function subs_todo = get_subs_list(indices)  

% SUBJECT EXCLUSION
% ------------------------ %
exclude_list = 31; % Due to brain abnormality

% warning('EXCLUDING THE 4 SBJS that have only 7 runs instead of 8')
% exclude_list = [exclude_list, 3, 13, 24, 33]; % due to only 7 runs
% 
% warning('EXCLUDING SBJ 1 that has 4 images less')
% exclude_list = [exclude_list, 1];
% ------------------------ %

disp('Calling subs_todo = get_subs_list(); will list all subjects for current analysis.')

global all_paths;
if isempty(all_paths)
    global_all_paths;
end
% Add SPM
addpath(all_paths.spm12)

% Get subs folders 
sub_dirs = cellstr(spm_select('FPList', all_paths.data_base_dir, 'dir', '^sub-[0-9]'));

% Parse subject number
subs_todo = zeros(size(sub_dirs'));
for s_ind = 1:length(sub_dirs)
    curr_sub = str2double(sub_dirs{s_ind}(end-1:end));
    if any(curr_sub == exclude_list)
        remove = 1;
        continue; 
    else
        subs_todo(s_ind) = curr_sub;
    end
end

% Remove all zeros from list
if remove == 1
    subs_todo = nonzeros(subs_todo)';
end

% Select subset of subjects
if exist('indices', 'var') && ~isempty(indices)
    subs_todo = subs_todo(indices);
end

fprintf('\nFound folders for n = %i subjects.\n', length(subs_todo))

end 