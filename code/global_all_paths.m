% WARNING! IF YOU CHANGE SOMETHING HERE, MAKE SURE TO COPY THE CHANGES TO
% ALL OTHER global_all_paths.m (if applicable)

function global_all_paths()
    % if global_all_paths_private.m is already in the matlab path
    if isempty(which('global_all_paths_private'))
        % check the private function is somewhere down the path
        p = fileparts(which(mfilename)); % get position of this file
        while 1
            lastp = p;
            p = fileparts(lastp); % go one level up
            if strcmp(p, lastp) % check if that was successful
                error('global_all_paths failed because could not find global_all_paths_private.m down to root, please check why')
            end
                
            % get position of this file
            if exist(fullfile(p, 'global_all_paths_private.m'), 'file')
                % add the directory above to the path
                display(['Adding global_all_paths_private.m from ' p])
                addpath(p)
                break
            end
        end
    end

    % call global_all_paths_private to really set the global path
    global_all_paths_private()
end