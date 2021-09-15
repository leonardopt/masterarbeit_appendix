% Start level1 FIR1 with condor
%
% to use condor:
% addpath('/analysis/share/condor_parallel/matlab_condor_parallel')
% substodo = 1:56
% % substodo([1, 3, 4, 5]) = [];
% create_parallel_kai('condor_level1_FIR1', substodo, '/home/kai/projects/aeffect_218/Analysis', 'condor_level1_FIR1')
%
% or manually (without logs)
% for job_nr = substodo
%   condor_level1_FIR1(job_nr)
% end

function condor_level1_FIR1_test(job_nr)
%% set path to your m-files

disp(['Starting script for job nr ' num2str(job_nr)])

global all_paths
%if isempty(all_paths)
    global_all_paths()
%end

% add current path
curr_path = pwd;
addpath(genpath(curr_path))  % adds path only for this session
% add SPM path
addpath(genpath(all_paths.spm12))

% !!! check on nathan, that path is set correct

%% copy your data from the server to the local machine

% if you have very many file accesses, please copy the file before running
% your script to the nathan, and copy the results back at the end
% (increases your speed as well as the speed for all others)
% 
% DONT FORGET TO CHANGE CHANGE THE PATH IN YOUR SCRIPT OR FUNCTION
% (e.g. by passing it as an argment, when you call it)

%% run your script

    % if you copy files, make sure to surround your script with a try/catch
    % loop, so that you can still copy the resulting files back
    try
        %%
        disp(['Running script for job nr ' num2str(job_nr)])
        
        subj = sprintf('%.2i',job_nr);
        
        inputs = {
            % BIDS root directory
            fullfile(all_paths.data_base_dir,'derivatives','preproc-1');
            % Output directory
            fullfile(all_paths.data_base_dir,'derivatives','level1-FIR1');
            % Participant or group level analysis
            'participant';
            % Flag to activate to run analysis on given participants, leave out for
            % all participants
            '--participant_label';
            % Labels of participants to include, leave out unless using
            % --participant_label flag
            subj;
            % Include for built-in BIDS validator, doesn't work properly, better to
            % use standalone version of bids-validator
            %'--bids-validator';
            % Keep original files
            %'--keep';
            % Indicate use of custom analysis pipeline
            '--config';
            % Path to pipeline script, default: pipeline_participant.m
            fullfile(all_paths.svn,'Analysis','level1','level1_FIR1.m')
            } 
        
        spm_BIDS_App
        disp(['Script for job nr ' num2str(job_nr) ' - ok'])

    catch my_error
        % display error and stack, and continue with the rest of the script
        display_error_info(my_error)
        disp(['Script for job nr ' num2str(job_nr) ' aborted'])
    end

%% copy you data back and clean up

% in case you are using a local directory, copy all data back to the net,
% and dont forget to delete the data on the remote machine afterwards

disp(['Finished script for job nr ' num2str(job_nr)])

end

function display_error_info(ME)
    ME
    for curr_field = fieldnames(ME)'
        disp([curr_field{1} ':'])
        display(ME.(curr_field{1}))
        disp(' ')
    end
    for stack_idx = 1:size(ME.stack, 1)
        display(ME.stack(stack_idx))
        disp(' ')
    end
end