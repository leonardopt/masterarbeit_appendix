%% Bonferroni-Holm corrections
%
% [H, nullrejections, info, belowchancerejects] = compute_bonf_holm_cor(p_values,
% alpha, belowchance, means, chancelevel)
%
% Bonferroni-Holm: Loop until ascending p-value exceeds
% alpha/(m-H0 rejected) count. Only accepts vectors. Expects uncorrected
% alpha (e.g. 0.05), not corrected values (but can of course be altered).
% If "belowchance" exists and is 0/false, below-chance accuracies will not
% be detected and added as H(ind) = -1; default is 1.
%
% H = vector or matrix with 1 for alternative hypothesis or 0 for null
% hypothesis
% nullrejections = total # of rejected null hypotheses for that vector
% info = short explanation
%
% Ingmar 09-09-2020

function [H, nullrejections, info, belowchancerejects] = compute_bonf_holm_cor(p_values, alpha, ...
    belowchance, means, chancelevel)

assert(length(size(p_values)) > 1, ...
    'compute_bonf_hom_cor() only accepts vectors of P-values');

info = [ 'Bonferroni-Holm correction incrementally decreases Type II error', ...
    ' after each rejected null hypothesis. No singular alpha or CI provided.'];

% Values
[pval_sort, sort_ind]   = sort(p_values); % sort
H                       = zeros(size(pval_sort)); % init H
m_tests                 = length(p_values); % init counter
nullrejections          = []; % init

% Above-chance sorting
for p_ind = 1:length(p_values)
    if pval_sort(p_ind) < (alpha/m_tests) % above-chance
        H(sort_ind(p_ind))  = 1; % index accordings to sort order
        m_tests             = m_tests - 1; % update m counter
        dispv(2, '    p = %.8f < a = %.8f', pval_sort(p_ind), alpha/m_tests);
    elseif pval_sort(p_ind) > (alpha/m_tests)
        nullrejections = p_ind - 1; % total # of null rejections per vector
        dispv(2, '    Holm-Bonf: procedure stopped at p = %.5f > a = %.5f, %i H0 rejected.', ...
            pval_sort(p_ind), alpha/m_tests, nullrejections);
        break; % stopping rule
    end
end % n_tests

if isempty(nullrejections) && p_ind == m_tests % stopping rule never reached
    nullrejections = p_ind;
    dispv(1, '    Holm-Bonf: all %i H0 rejected.', m_tests);
end

% Below-chance sorting (additional test; default = 1)
if exist('belowchance', 'var') && belowchance == 0
    warning('Flag "belowchance" == 0, quitting MCC correction')
    return; % do not continue significant below-chance accuracy detection
end

assert(exist('means', 'var') && ~isempty(means) && ...
    exist('chancelevel', 'var') && ~isempty(chancelevel),...
    'Input means or chancelevel is missing or empty.');

b_tests     = length(p_values); % init counter
bval_sort   = flip(pval_sort);
bval_ind    = flip(sort_ind);

for b_ind = 1:b_tests
    if bval_sort(b_ind) > 1-(alpha/b_tests)
        H(bval_ind(b_ind))  = -1; % index accordings to sort order
        b_tests             = b_tests - 1; % update b counter
        dispv(2, '    p = %.8f > 1-a = %.8f', pval_sort(p_ind), alpha/b_tests);
    else
        belowchancerejects = b_ind - 1;
        dispv(2, '    Holm-Bonf: procedure stopped at p = %.5f < 1-a = %.5f, %i H0 rejected.', ...
            bval_sort(b_ind), alpha/b_tests, belowchancerejects);
        break; % stopping rule
    end
end % n_tests

isnull = squeeze(round(p_values, 4) == 1) & squeeze(round(means, 4) == ...
    round(chancelevel, 4)); % equal to chance
H(isnull == 1) = 0;
dispv(2, '    Holm-Bonf: corrected %i false below-chance testing outcomes to zero', sum(isnull));

end % func