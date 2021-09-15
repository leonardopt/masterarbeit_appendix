% Function to find minimum offset to add to the data, in order to avoid
% spm_global unwanted masking of voxels 
%
% See spm_global for details.
% 
% Specifically, spm_global checks that 
% ALL VALUES IN ALL VOXELS that should not be masked out MUST BE > GLOBAL MEAN / 8
%
% Here the excerpt from the spm_global header
%
% ---
% % spm_global returns the mean counts integrated over all the slices from
% the volume.
%
% The mean is estimated after discounting voxels outside the object using
% a criteria of greater than > (global mean)/8.
% --
%
% Kai and Leonardo, some time in 2019/20

function target_offset = find_offset_for_SPM_global_voxelvalues(dat)

alldat = vertcat(dat(:).X);

% some checks
n_nans = isnan(alldat(:));
if n_nans == numel(alldat)
    error('find_offset_for_SPM_global_voxelvalues:all_nan', ...
        'All data is nan, cannot make it SPM compatible. consider replacing it with all 0s if its necessary.')
elseif n_nans > 0
    warning('find_offset_for_SPM_global_voxelvalues:some_nan', ...
        '%i/%i values are nan. Trying to make the data spm compatible. Not sure if it works.', ...
        n_nans, numel(alldat))
end

fglobalmean = @(x)nanmean(alldat(:) + x)/8; % x is intercept 
fminvals  = @(x)nanmin(alldat(:) + x);

% Create vector with potential intercept values (ranging from 0 to n)
xintercepts = 0:1:10;

% Get vector of mean/min values given each potential intercept  
globalmeanvals     = fglobalmean(xintercepts);

i = 1;
stopcounter = 10000;
counter = 0;
while ~( fminvals(xintercepts(i)) > globalmeanvals(i) )
    i = i + 1;
    if i > numel(xintercepts)
        xintercepts = xintercepts + xintercepts(end);
        i = 1; % restart counter 
    end     

    % assert that no compontent is nan
    assert(~isnan(fminvals(xintercepts(i))), 'xintercept is nan but should not be');
    assert(~isnan(globalmeanvals(i)), 'globalmeanvals is nan but should not be');

    % check for number of iterations
    counter = counter + 1;
    assert(stopcounter > counter, 'find_offset_for_SPM_global_voxelvalues:while_loop_limit_exceeded', ...
        'Could not find a value in reasonable time. check if there is a problem in the data.');
end 

target_offset = xintercepts(i);

end 