% function [nii_dat, nii_nam] = create_nifti_multiple_models(dats_to_write, fname)
%
% Function that creates a nifti image
%
% IN
% dats_to_write: data that will be written in the image
% fname: name of the output file
%
% OUT
% nii_dat: voxel data
% nii_nam: voxel names
%
% NB: each full model is written in the first two dimensions of the image.
% Different full models are written in the third dimension.
% e.g. if dats_to_write is data of five full models, each one
% containing 9 voxels and 200 scans, the nifti image will have
% dimensionality 3 x 3 x 5 x 200
%
% NB2: if only one full model is present, the scripts adds a fake second
% model with all nans, to ensure the presence of a third dimension
% (otherwise SPM might trip up)
%
% Kai/Leonardo, 2019

function [nii_dat, nii_nam] = create_nifti_multiple_models(dats_to_write, fname, add_fullmodellayer_to_voxelname)

%% defaults

if ~exist('add_fullmodellayer_to_voxelname', 'var')
    add_fullmodellayer_to_voxelname = 1; % fallback to old default behaviour
end

%% File name

% if file name does not have .nii extension, add it
if  ~contains(fname, '.nii')
    fname = [fname '.nii'];
end

%% If dimensionality of dats to write is one, add a second empty model to ensure 3 dimensionality

if length(dats_to_write) == 1
    dats_to_write(2).X = nan(size(dats_to_write(1).X));
    dats_to_write(2).voxelnames = cell(size(dats_to_write(1).voxelnames));
    if isfield(dats_to_write, 'weights')
        dats_to_write(2).weights = zeros(size(dats_to_write(1).weights));
    end
end

if add_fullmodellayer_to_voxelname
    for d_ind = 1:length(dats_to_write)
        dats_to_write(d_ind).voxelnames = strcat(dats_to_write(d_ind).voxelnames, sprintf('_fullmodel%i', d_ind));
    end
end

%% Define the four dimensions

nregressors = size(dats_to_write(1).X, 2);

dim1 = ceil(sqrt(nregressors)); % ceiling of sqrt of n_vox
dim2 = dim1;
dim3 = length(dats_to_write); % n of models
dim4 = size([dats_to_write.X], 1);

% Define target shape of the nifti image
nii_dat = nan(dim1, dim2, dim3, dim4); % pre-allocate

%% Adjust dimensionality by adding nans if necessary

nanstoadd = dim1*dim2 - nregressors;
if nanstoadd ~= 0
    for d_ind = 1:length(dats_to_write)
        dats_to_write(d_ind).X          = [dats_to_write(d_ind).X nan(size(dats_to_write(d_ind).X, 1), nanstoadd)];
        dats_to_write(d_ind).voxelnames = [dats_to_write(d_ind).voxelnames cell(size(dats_to_write(d_ind).voxelnames,1), nanstoadd)];
    end

    % Check equal voxel sizes in 3rd dim
    err_str = ['Creating nifti models is (currently) incompatible with different voxel sizes per layer. ', ...
        'Voxel count needs to be consistent across the 3rd dimension.'];
    assert((size(nii_dat, 1) * size(nii_dat, 2) * size(nii_dat, 3)) == size([dats_to_write.X], 2), err_str);
end

%% Get data from dat file

alldat     = [dats_to_write.X];
voxelnames = [dats_to_write.voxelnames];
voxelnames = repmat(voxelnames, dim4, 1);

%% Reshape data in 4D

nii_dat = reshape(alldat', size(nii_dat));
nii_nam = reshape(voxelnames', size(nii_dat));

%% Create image

ni           = nifti;
% header info
ni.mat       = eye(4);
ni.dat       = file_array(fname, size(nii_dat), [0 spm_platform('bigend')], 0, 1, 0); % SPM
ni.descrip   = fname; % CHANGE THIS
ni.dat.dtype = 'FLOAT32';
create(ni);

%% Write image

try
    for d_ind = 1:size(ni.dat,4)
        ni.dat(:,:,:,d_ind) = nii_dat(:,:,:,d_ind);
        spm_get_space([ni.dat.fname ',' num2str(d_ind)], ni.mat); % SPM
    end
catch writingerror
    disp(['    Error: ' writingerror.message]);
    error('Writing nifti images failed (possibly due to writing permission denied)...');
end

end

% % If optional input is a struct, turn it into a string
% if strcmp(class(description), 'struct')
%     if length(fieldnames(description)) == 2 && sum(ismember(fieldnames(description), 'name')) == 1 && sum(ismember(fieldnames(description), 'voxels'))
%                 description_string = sprintf('[condition: voxel num]\n');
%                 for i = 1:length(description)
%                     description_string = sprintf([description_string description(i).name ': ' num2str(description(i).voxels) '\n']);
%                 end
%                 description_string(end) = ''; % get rid of last linebreak
%                 description = description_string;
%
%     else
%         error('Input must be a struct with two fields called ''name'' and ''voxels''.');
%     end
% end

% 'descrip' field of a header cannot be longer than 80 characters
% if 'description' variable is longer than 80 char, write it in an external
% file
% if length(description) > 80
%     info_fname = [fname(1:end-4) '_info']; % remove .nii extension
%     save(info_fname, 'description');
%     description = 'Description in _info file';
%     fprintf('Description length greater than 80 characters. Saving info file...\nInfo file directory: %s...\n', info_fname);
% end
