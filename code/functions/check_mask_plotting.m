% [fig_title, maskfile] = check_mask_plotting(result, fig_title)
%
% Separated from plot_decoding_results. Checks mask info before putting in
% (sub)plot. Gets maskfile of first + test subject.
%
% Kai, September 2021

function [fig_title, maskfile] = check_mask_plotting(result, fig_title)

if ~exist('fig_title', 'var')
    try
        fig_title = [result.mask_name ' s' num2str(result.subs_todo(1))];
    catch
        fig_title = 'Current mask';
    end
end

% Get mask
disp(['Loading ' result.full_datafile{1} ' to get maskfile to plot'])
c1 = load(result.full_datafile{1});
try
    maskfile = c1.ACC_mask_data.maskfile;
catch
    maskfile = c1.ACC_mask_data{1}.maskfile;
end
[~, fn, ~] = fileparts(maskfile);
if ~strcmp(result.mask_name, fn)
    disp('    Maskfile in loaded data does not fit to result.mask_name. This can be an old bug before recomputing. Checking if loaded .info has the expected maskname, and if so, changing the filename')
    if strfind(c1.ACC_mask_data.info, result.mask_name)
        disp('    Found expected maskname in string, replacing filename to load')
        maskfile = strrep(maskfile, fn, result.mask_name);
    else
        error('Maskfile in loaded data and in result.mask_name differ, please check')
    end
end
disp(['    Found maskfile for plotting: ' maskfile])

% Check mask
checksbj_ind = randi([2, length(result.full_datafile)]);
disp(['  Loading data from random subject ' result.full_datafile{checksbj_ind} ' to get maskfile to check'])
c2 = load(result.full_datafile{checksbj_ind});
try
    checkmaskfile = c2.ACC_mask_data.maskfile;
catch
    checkmaskfile = c2.ACC_mask_data{1}.maskfile;
end
[~, fn, ~] = fileparts(checkmaskfile);
if ~strcmp(result.mask_name, fn)
    disp('    Checkmaskfile in loaded data does not fit to result.mask_name. This can be an old bug before recomputing. Checking if loaded .info has the expected maskname, and if so, changing the filename')
    if strfind(c2.ACC_mask_data.info, result.mask_name)
        disp('    Found expected maskname in string, replacing filename to load')
        checkmaskfile = strrep(checkmaskfile, fn, result.mask_name);
    else
        error('Maskfile in loaded data and in result.mask_name differ, please check')
    end
end
disp(['    Found maskfile for checking: ' checkmaskfile])

% Read and sum mask
maskdata  = spm_read_vols(spm_vol(maskfile));
checkdata = spm_read_vols(spm_vol(checkmaskfile));
summask = sum(maskdata(:)>0);
sumcheck = sum(checkdata(:)>0);

% Check sum
if ~isequal(size(maskfile), size(sumcheck))
    warning('plot_decoding_results:subject_masks_differ_dimension', ...
        'Subject masks differ in dimension, ok for decoding in subjects space')
    fig_title = [fig_title '(masksize differs)'];
else
    sumdiff = sum(maskdata(:)>0~=checkdata(:)>0);
    if sumdiff > 0
        warning('plot_decoding_results:subject_masks_differ', ...
            'Subject masks have same dimension but data differ: sum(mask(:)>0)=%i vxls, sum(check(:)>0)=%i vxls, sumdiff (mask(:)>0~=check(:)>0): %i vxls', summask, sumcheck, sumdiff)
    else
        disp('    Both masks agree.')
    end
    fig_title = [fig_title '(maskcontent differs)'];
end

end