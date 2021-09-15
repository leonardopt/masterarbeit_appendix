% function dat = add_voxels_with_mixed_effects(dat, n_mixedvoxels, weight)
%
% Adds voxels with mixed (weighted) responses to task/rest structure found
% in the empirical data. Function as noise voxels. Can be tweaked using
% simulation parameters found in the setup respective to the simulation
% name.

function data = add_voxels_with_mixed_effects(data, n_mixedvoxels, weight)

if ~exist('weight', 'var')
    weight = 1;
end

%% Add voxels with mixed effects

for r = 1:length(data)
    currX = data(r).X;
    for i = 1:n_mixedvoxels
        voxel       = currX * randn(size(currX, 2), 1) * weight;
        data(r).X   = [data(r).X voxel]; % add voxel to design mat
        data(r).voxelnames{end+1} = sprintf('mixed_effects_voxel_%i', i)  ;% add voxel name to voxel names field
    end
end

dispv(1, '    %i mixed response voxels added...', n_mixedvoxels);

end