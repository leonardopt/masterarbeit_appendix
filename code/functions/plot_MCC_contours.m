% sp_hdl = plot_MCC_contours(result, sp_hdl, customticks, customlabels)
%
% Plots all contours (outlines) of significance values of all possible
% multiple comparison corrections. As of 08-09-21, there are 6 plus
% uncorrected). These are:
%
% 0: uncorr
% 1: Bonferroni corrected
% 2: FDR
% 3: Sign permutation (clustersize) reconstructed from individual quandrant matrices
% 4: Sign permutation (clustermass) reconstructed from individual quandrant matrices
% 5: TFCE
% 6: TFCE reconstructed from individual quandrant matrices
%
% As background, it puts a greyscale version of the accuracy matrix.
%
% Reliant on results struct from afx_grouplevel_decoding_stats and subplot
% handle but could easily be modified as standalone figure too.
%
% Ingmar 08-09-21

function sp_hdl = plot_MCC_contours(result, sp_hdl, customticks, customlabels)
%% Collect subplot MCC stats

% Uncorr
MCC_stats = {}; % init
MCC_stats{1} = result.uncorr.H;

% Bonferroni corrected
MCC_stats{2} = result.bonf.H;

% FDR
MCC_stats{3} = result.fdr.P;

% Sign permutation (clustersize) reconstructed from individual quandrant matrices
if isfield(result.signperm, 'dq') && length(result.signperm.dq) >= 6 && isfield(result.signperm.dq(6), 'H_size')
    MCC_stats{4} = result.signperm.dq(6).H_size;
elseif isfield(result.signperm, 'H_size')
    MCC_stats{4} = result.signperm.H_size;
    % old format
elseif isfield(result.signperm, 'dq') && length(result.signperm.dq) >= 6 && isfield(result.signperm.dq(6), 'H')
    warning('old signperm format, is onlz .H but should now be .H_size')
    MCC_stats{4} = result.signperm.dq(6).H;
elseif isfield(result.signperm, 'H')
    warning('old signperm format, is onlz .H but should now be .H_size')
    MCC_stats{4} = result.signperm.H;
else
    warning('Could not find data for sign permutation (clustersize)')
    MCC_stats{4} = []; %zeros(size(result.mean)); % empty array
end

%  Sign permutation (clustermass) reconstructed from individual quandrant matrices
try
    MCC_stats{5} = result.signperm.dq(6).H_mass;
catch
    if isfield(result.signperm, 'H_mass')
        MCC_stats{5} = result.signperm.H_mass;
    elseif isfield(result.signperm, 'dq') && length(result.signperm.dq) >= 6 && isfield(result.signperm.dq(6), 'H_size') && ... % H_size correct here, only to check for format (the mass data is in H)
            isfield(result.signperm.dq(6), 'H') && isfield(result.signperm.dq(6), 'name') && strfind(result.signperm.dq(6).name, '(mass)')
        warning('sign permutation (clustermass): Using temporary signperm.H because (mass) occurs in name')
        MCC_stats{5} = result.signperm.dq(6).H;
    else
        warning('Could not find data for sign permutation (clustermass)')
        MCC_stats{5} = []; % zeros(size(result.mean)); % empty array
    end
end

% TFCE
try
    P_neg = result.tfce.dq(1).pcorr_neg2;
    P_neg(P_neg==1) = .5;
    MCC_stats{6} = P_neg;
catch
    warning('Could not find data for tfce permutation')
    MCC_stats{5} = []; %zeros(size(result.mean)); % empty array
end

% TFCE reconstructed from individual quandrant matrices
try
    MCC_stats{7} = result.tfce.dq(6).H;
catch
    warning('Could not find data for tfce permutation with dq reconstruction')
    MCC_stats{7} = []; %zeros(size(result.mean)); % empty array
end

%% Plotting

% Correction for matlab axes
offset = 0.5;

% Show mean in grayscale
%figure('name', 'Multiple comparison corrections (MCC) contour lines for significance
hold on; shading('flat'); imagesc(sp_hdl, result.mean(1:43, 1:43));
colormap(sp_hdl, 'gray');

% Colors
[~, ~, plasma, viridis] = imported_colormaps;
neg_ind = 5:30:5+(30*7);
neg = viridis(neg_ind, :);
pos_ind = 240:-15:240-(15*7);
pos = plasma(pos_ind, :);

% Show pos contours
for s_ind = 1:length(MCC_stats)
    if ~isempty(MCC_stats{s_ind})
        Z = MCC_stats{s_ind}(1:43, 1:43);
        [Y,X] = meshgrid(1:size(Z, 1), 1:size(Z,2)); % Decrease line overlap
        contour(X' + s_ind * -0.1, Y' + s_ind * -0.1, Z, [0.99 1], 'LineColor', pos(s_ind, :), 'LineWidth', 1);
    end
end

% Show neg contours
for s_ind = 1:length(MCC_stats)
    if ~isempty(MCC_stats{s_ind})
        Z = MCC_stats{s_ind}(1:43, 1:43);
        [Y,X] = meshgrid(1:size(Z, 1), 1:size(Z,2)); % Decrease line overlap
        contour(X' + s_ind * -0.1, Y' + s_ind * -0.1, Z, [-1 -0.99] , 'LineColor', neg(s_ind, :), 'LineWidth', 1);
    end
end
hold off;

% Names for methods present (no empty stats)
names = {'Uncorrected', 'Bonferroni', 'FDR 0.001', 'Sign perm. (size)', ...
    'Sign perm. (mass)', 'TFCE', 'TFCE quadr.'};
names_ind = ~cellfun(@isempty, MCC_stats);
names = names(names_ind);
allnames = cat(2, names, names);
l_hdl = legend(allnames, 'Location', 'eastoutside', 'AutoUpdate', 'off', 'FontSize', 8);
try l_hdl.ItemTokenSize = [10 10 10]; end % reduce legend line length

% Rest phase lines
cline = [0 0 0]; % [0.6 0.6 0.6];
xline(27 + offset, 'Color', cline, 'LineWidth', 1); % 'LineStyle', '--');
yline(27 + offset, 'Color', cline, 'LineWidth', 1); % 'LineStyle', '--');

% Alter color limits
if isfield(result, 'extend_caxis') && ~isempty(result.extend_caxis)
    caxis(result.extend_caxis(1, :));
end

% Add task lines
xline(27 + offset, 'w-', 'LineWidth', 1);
yline(27 + offset, 'w-', 'LineWidth', 1);
hold off;

% add xlabel for accuracies
xlabel(sp_hdl, 'Decoding accuracy (%)');

% Get axes ticks and labels
title('X. MCC contours', 'Interpreter', 'none');
axis equal; axis tight;
xlabel('Testing time-points'); ylabel('Training time-points');
set(gca, 'XTickLabelRotation', 60); set(gca,'FontSize',8);
yticks(customticks + offset); xticks(customticks + offset);
yticklabels(customlabels); xticklabels(customlabels);
set(gca, 'TickDir', 'out');

end