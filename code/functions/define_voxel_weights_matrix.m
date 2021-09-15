% function voxels = define_voxel_weights_matrix(pw, condweights, taskweights)
%
% Convolves generated voxels with certain weight parameters set prior:
% pw (prior weighting), weights separate for each phase (cue, task, rest),
% and includes separate weights for the three tasks.
%
% Leonardo, 2019

function voxels = define_voxel_weights_matrix(pw, phaseweights, taskweights, normalisemat)
%% Preferred condition

%              Aeh | Kon | Sum
voxels =    [[ 1*pw   1     1    ]; % voxel 1 - cue
            [  1      1*pw  1    ]; % voxel 2 - cue
            [  1      1     1*pw ]; % voxel 3 - cue
            [  1*pw   1     1    ]; % voxel 4 - trial
            [  1      1*pw  1    ]; % voxel 5 - trial
            [  1      1     1*pw ]; % voxel 6 - trial
            [  1*pw   1     1    ]; % voxel 7 - rest
            [  1      1*pw  1    ]; % voxel 8 - rest
            [  1      1     1*pw ]; % voxel 9 - rest
            ];

%% Apply additional weights

% 1) By condition
voxels(1:3, :) = voxels(1:3, :) * phaseweights(1); % cueweight
voxels(4:6, :) = voxels(4:6, :) * phaseweights(2); % trialweight
voxels(7:9, :) = voxels(7:9, :) * phaseweights(3); % restingweight

% 2) By task
voxels(:, 1) = voxels(:, 1) * taskweights(1); % aehweight;
voxels(:, 2) = voxels(:, 2) * taskweights(2); % konweight;
voxels(:, 3) = voxels(:, 3) * taskweights(3); % sumweight;

%% Normalise rows if selected (warning: will make all data very similar)

if exist('normalisemat', 'var') && normalisemat == 1
    for i = 1:size(voxels,1)
        voxels(i, :) =  voxels(i, :)/sum(voxels(i, :));
    end
end

end