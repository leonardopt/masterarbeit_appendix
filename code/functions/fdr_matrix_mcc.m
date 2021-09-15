% result = fdr_matrix_mcc(data, chancelvl, base_alpha, FDR_alpha, sz, means, nodisp)
%
% Corrects 2D matrix of data (eg. timeresolved accuracy matrix) for
% multiple comparisons using the Benjamini and Hochberg 1995 false discovery
% rate algorithm with the potential assumption of dependent data.
%
% Needs fdr_bh (included) but a possible fallback is the matlab option 
% mafdr(old code included below)
%
% Note: CI95 not included since they would only be generated for the
% significant values.
%
% Ingmar, May 2021

function result = fdr_matrix_mcc(result, chancelvl, base_alpha, FDR_alpha, sz, means, nodisp)
%% Checks

% For below-chance hypothesis check
if ~isfield(result, 'mean') 
    if logical(exist('mean','var')) && ~isempty(means)
        result.mean = means;
    else
         error('Mean original data values not supplied with results.mean or as last variable, please add them!');
    end    
end
    
% Compatibility
if ~isstruct(result)
    notstruct = result; result = struct;
    result.uncorr.P = notstruct; % must be P-values
elseif ~isfield(result, 'uncorr') || ~isfield(result.uncorr, 'P')
    try
        result.uncorr.P = result.P;
    catch
        error('No P-values found at assumed locations! Provide matrix of P-values as result.uncorr.P');
    end
end
assert(numel(size(result.uncorr.P)) == 2, 'Matrix is not 2D!');
if ~exist('size', 'var') || isempty(sz)
    sz = size(result.uncorr.P, 1);
end

% Chance lvl
if ~exist('chancelvl', 'var') || isempty(chancelvl)
    chancelvl = 100/3;
end

% Alpha lvl
if ~exist('base_alpha', 'var') || isempty(base_alpha)
    base_alpha = 0.05;
end

% FDR
if ~exist('FDR_alpha', 'var') || isempty(FDR_alpha)
    FDR_alpha = 0.001;
end

% Disp output
if ~exist('nodisp', 'var') || isempty(nodisp)
    nodisp = 'yes';
    fprintf('\n    Matrix size for FDR is set as %i x %i\n', sz, sz);
    fprintf('\n    FDR %f correction of two-sided alpha = 0.05\n', FDR_alpha);
elseif nodisp == 1
    nodisp = 'no';
end

%% B&H algorithm, could also yield CI but only for sig. vals

[result.fdr.H, ~, ~, result.fdr.P] = fdr_bh(result.uncorr.P(:), FDR_alpha, 'pdep', nodisp);
result.fdr.P = reshape(result.fdr.P, size(result.uncorr.P));
result.fdr.H = double(result.fdr.P < base_alpha); % pos

% below chance
[~, ~, ~, result.fdr.P_below] = fdr_bh(1-result.uncorr.P(:), 0.001, 'pdep', nodisp);
result.fdr.P_below = reshape(result.fdr.P_below, size(result.uncorr.P));
result.fdr.H(result.fdr.P_below < base_alpha) = -1;
isnull = round(result.fdr.P, 4)==1 & (round(result.mean, 4)==round(chancelvl, 4) == 1); % equal to chance
result.fdr.H(isnull == 1) = 0;

end


%% OLD FDR CODE
%
%     % Normal above chance
%     [result.fdr.P, result.fdr.error, result.fdr.apriori] = mafdr(result.uncorr.P(:));
%     result.fdr.P = reshape(result.fdr.P, size(result.uncorr.P));
%     result.fdr.error = reshape(result.fdr.error, size(result.uncorr.P));
%     result.fdr.H = double(result.fdr.P < base_alpha/2); % pos
%
%     % Normal below chance
%     [result.fdr.P_below, result.fdr.error_below, result.fdr.apriori_below] = mafdr(1-result.uncorr.P(:));
%     result.fdr.P_below = reshape(result.fdr.P_below, size(result.uncorr.P));
%     result.fdr.H(result.fdr.P_below < base_alpha/2) = -1;
%     isnull = round(result.fdr.P, 4)==1 & (round(result.mean, 4)==round(chancelvl, 4) == 1); % equal to chance
%     result.fdr.H(isnull == 1) = 0;

%
%     for c_ind = 1:length(comp_results)
%         assert(isfield(comp_results{c_ind}, 'ACC_elements') && ~isempty(comp_results{c_ind}), ...
%             'ACC_elements missing: include 4-D array of accuracy matrices over subjects (dim 3) and comparisons (dim 4)');
%         curr_ACC_elements(:, :, :, c_ind) = comp_results{c_ind}.ACC_elements;
%     end
%