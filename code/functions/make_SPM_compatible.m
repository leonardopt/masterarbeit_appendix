% function [all_data, make_SPM_compatible_info] = make_SPM_compatible(all_data)
%
% To be applied BEFORE spm_global to avoid masking in spm_global.
%
% Made by Kai and Leonardo, 2019. Edited for compatibility with multiple
% nifti layers, Ingmar Jan 2021

function [all_data, make_SPM_compatible_info] = make_SPM_compatible(all_data)

% Get minimum value for the current subject
for d_ind = 1:length(all_data)    
    % takes minimum offset that should be added to the data to satisfy the condition: min(dat) > mean(dat)/8
    target_offset(d_ind) = find_offset_for_SPM_global_voxelvalues(all_data{d_ind});  %#ok<AGROW>
end

% Get maximum of all minima to equalize offsets across layers
max_target_offset = max(target_offset);


% Determine class
cl = class(all_data{1});

switch cl
    case 'struct'
                
        for d_ind = 1:length(all_data)

            curr_layer = all_data{d_ind}; % unpack
            
            % layer-wise
            for l_ind = 1:length(curr_layer)
                curr_layer(l_ind).X = curr_layer(l_ind).X + max_target_offset;
            end
            all_data{d_ind} = curr_layer; % put back
        end
        
    case 'double'
        
        for d_ind = 1:length(all_data)
            
            % Add offset
            curr_layer = all_data{d_ind};
            curr_layer = curr_layer + max_target_offset;
            all_data{d_ind} = curr_layer; % put back
        end
    otherwise
        error('Variable "all_data": check class')
end

make_SPM_compatible_info.offset = max_target_offset;

dispv(1, '    Adding offset of %g to all data to make it SPM compatible (avoids masking in spm_global).', max_target_offset)