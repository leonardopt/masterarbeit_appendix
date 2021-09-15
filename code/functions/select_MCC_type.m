% [subplot_stat, subplot_titles, MCC_type, targetfile] = select_MCC_type(curr_MCC, result, targetfile)
%
% Selects desired MCC statistic from result struct.
%
% Separated from plot_decoding_results.m
%
% Ingmar and Kai 10-09-21

function [subplot_stat, subplot_titles, MCC_type, targetfile] = select_MCC_type(curr_MCC, result, targetfile)

% init
subplot_stat = {};
subplot_titles = {};
MCC_type = '';

if ~exist('curr_MCC', 'var') || ~curr_MCC % 0
    subplot_stat    = {result.mean, result.uncorr.CItwos_diff, result.uncorr.H};
    subplot_titles  = {'A. Mean accuracies', 'B. CI95-width (uncorr.)', ...
        'C. Significance (uncorr.)'};
    % Save maps, add mask name and comparison condition to output file
    MCC_type = ' uncorr.';
    targetfile = [targetfile '_uncorr'];
    
elseif curr_MCC == 1 % bonf
    subplot_stat = {result.mean, result.bonf.CItwosided_diff, result.bonf.H};
    subplot_titles  = {'A. Mean accuracies', 'B. CI95-width (Bonf.)', ...
        'C. Significance (Bonf.)'};
    MCC_type = ' Bonf. corr.';
    targetfile = [targetfile '_bonf'];
    
elseif curr_MCC == 2 % fdr
    subplot_stat = {result.mean, result.uncorr.CItwos_diff, result.fdr.H, result.fdr.P};
    subplot_titles  = {'A. Mean accuracies', 'B. CI95-width (uncorr.)', ...
        'C. Significance (FDR 0.001)', 'X. P_FDR'};
    MCC_type = ' FDR';
    targetfile = [targetfile '_FDR'];
    
elseif curr_MCC == 3 % sign permutation (clustersize)
    % size and mass format
    if isfield(result.signperm, 'dq') && length(result.signperm.dq) >= 6 && isfield(result.signperm.dq(6), 'H_size')
        subplot_stat = {result.mean, result.uncorr.CItwos_diff, result.signperm.dq(6).H_size};
    elseif isfield(result.signperm, 'H_size')
        subplot_stat = {result.mean, result.uncorr.CItwos_diff, result.signperm.H_size};
        % old format
    elseif isfield(result.signperm, 'dq') && length(result.signperm.dq) >= 6 && isfield(result.signperm.dq(6), 'H')
        warning('old signperm format, is onlz .H but should now be .H_size')
        subplot_stat = {result.mean, result.uncorr.CItwos_diff, result.signperm.dq(6).H};
    elseif isfield(result.signperm, 'H')
        warning('old signperm format, is onlz .H but should now be .H_size')
        subplot_stat = {result.mean, result.uncorr.CItwos_diff, result.signperm.H};
    else
        warning('Could not find data for sign permutation (clustersize), skipping')
        return
    end
    if isempty(subplot_stat{3})
        warning('sign permutation (clustersize): data was empty, skipping')
        return
    end
    subplot_titles  = {'A. Mean accuracies', 'B. CI95-width (uncorr.)', ...
        'C. Significance (Sign perm. size)'};
    MCC_type = ' Sign perm. clustersize';
    targetfile = [targetfile '_signperm_size'];
    
elseif curr_MCC == 4 % sign permutation (clustermass)
    try
        subplot_stat = {result.mean, result.uncorr.CItwos_diff, result.signperm.dq(6).H_mass};
    catch
        if isfield(result.signperm, 'H_mass')
            subplot_stat = {result.mean, result.uncorr.CItwos_diff, result.signperm.H_mass};
        elseif isfield(result.signperm, 'dq') && length(result.signperm.dq) >= 6 && isfield(result.signperm.dq(6), 'H_size') && ... % H_size correct here, only to check for format (the mass data is in H)
                isfield(result.signperm.dq(6), 'H') && isfield(result.signperm.dq(6), 'name') && strfind(result.signperm.dq(6).name, '(mass)')
            warning('sign permutation (clustermass): Using temporary signperm.H because (mass) occurs in name')
            subplot_stat = {result.mean, result.uncorr.CItwos_diff, result.signperm.dq(6).H};
        else
            warning('Could not find data for sign permutation (clustermass), skipping')
            return
        end
    end
    subplot_titles  = {'A. Mean accuracies', 'B. CI95-width (uncorr.)', ...
        'C. Significance (Sign perm. mass)'};
    MCC_type = ' Sign perm. clustermass';
    targetfile = [targetfile '_signperm_mass'];
    
elseif isfield(result, 'tfce') && curr_MCC == 5 % TFCE (full matrix, no reconsclctruction)
    P_pos = result.tfce.dq(1).pcorr_pos2;
    P_pos(P_pos == 1) = .5;
    P_neg = result.tfce.dq(1).pcorr_neg2;
    P_neg(P_neg==1) = .5;
    subplot_stat = {result.mean, P_pos, result.tfce.dq(1).H, P_neg};
    subplot_titles  = {'A. Mean accuracies', 'X. TFCE P_pos_', ...
        'C. Significance (TFCE, standard)', 'X. TFCE P_neg'};
    MCC_type = ' TFCE';
    targetfile = [targetfile '_tfce'];
    
elseif isfield(result, 'tfce') && curr_MCC == 6 % TFCE (reconstructed from submatrices)
    subplot_stat = {result.mean, result.uncorr.CItwos_diff, result.tfce.dq(6).H};
    subplot_titles  = {'A. Mean accuracies', 'B. CI95-width (uncorr.)', ...
        'C. Significance (TFCE, quadrantwise.)'};
    MCC_type = ' TFCE (submat recon)';
    targetfile = [targetfile '_tfcesubmatrecon'];
else
    warning('Unknown do_MCC number'); return
end

fprintf('\n   Multiple comparison method #%i "%s" data set for plotting\n', curr_MCC, subplot_titles{3});

end