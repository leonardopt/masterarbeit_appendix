% fighdl = plot_timecourse_CI_filled(result, customticks, ...
%    customlabels, do_MCC, do_diagonal, slice_vals, fighdl)
%
% Constructs plotting with timecourse and CI, color filled instead of bars,
% and saves it to specified folder. Also, it constructs slices of the training set with
% timecourse and CI, color filled, and saves it to specified folder. Uses
% two training timepoints in TR, e.g. Trial10 and Rest10 (or otherwise specified)
% instead of the full timecourse to show variability across test time.
%
% Merged: Ingmar 15/09/2020
% Made: Ingmar 02/06/2020

function fighdl = plot_timecourse_CI_filled(result, customticks, ...
    customlabels, do_MCC, do_diagonal, slice_vals, fighdl)

if ~exist('fighdl', 'var') || isempty(fighdl)
    fighdl = figure;
end

% Values
xaxis_values    = result.axislabels.testtime_xaxis_vals;
chancelevel     = result.chancelevel;
reststart       = 27;
offset          = 0.5;

% Labels
y_labels    = 'Accuracy (%) + CI95';
x_labels    = 'Testing time-points';


% Stats to show
if ~exist('do_MCC', 'var') || do_MCC ~= 1
    CItwosided = result.uncorr.CItwosided; % Uncorr CI for all methods except bonferroni
elseif do_MCC == 1
    CItwosided = result.bonf.CItwosided;
end
[subplot_stat, ~, ~, ~] = select_MCC_type(do_MCC, result, '');
H_vals = subplot_stat{3};

% Slice input
time_vals = slice_vals;

%% Constructing significance line

% COLORS
blue    = [53/255 119/255 233/255];
purple  = [105/255 55/255 168/255];
pinkish = [243/255 80/255 122/255];

switch do_diagonal
    case true
        %% DIAGONAL
        
        CI_lowbound     = diag(squeeze(CItwosided(1, :, :)));
        CI_highbound    = diag(squeeze(CItwosided(2, :, :)));
        
        if do_MCC == 1
            sig_vals = diag(H_vals);
        else
            sig_vals = diag(H_vals);
        end
        sig_pts = round(max(CI_highbound)+4) * ones(size(CI_highbound))';
        sig_pts(~sig_vals)      = NaN;
        sig_first               = find(~isnan(sig_pts), 1, 'first');
        sig_last                = find(~isnan(sig_pts), 1, 'last');
        dispv(1, '    Diagonal timecourse first and last significant values: %i, %i.', sig_first, sig_last);
        
        hold on; % Plot
        plot(xaxis_values, diag(result.mean), 'b-', 'LineWidth', 1);
        plot(sig_pts, 'k-', 'LineWidth', 1.3); % significance line
        
        % fill needs all edges defined (e.g. for polygons), so use fliplr
        % to give axis points in all "directions". not sure why flipud is
        % needed instead of fliplr... but okay?
        xaxis_fill = [xaxis_values, fliplr(xaxis_values)];
        yaxis_fill = [CI_highbound; flipud(CI_lowbound)];
        fill(xaxis_fill, yaxis_fill, blue, 'facecolor',  ...
            [53/255 119/255 233/255], 'edgecolor', [53/255 119/255 233/255], 'facealpha', 0.15);
        yline(chancelevel + offset + offset, 'k-'); % adds chance level line
        xline(reststart + offset + offset, 'k--'); % adds end of task phase
        hold off;
        
        % add legend to make plots more aligned
        legend('Diagonal', 'Location', 'southoutside', 'AutoUpdate', 'off');
        
        newtitle = 'D. Diagonal time-course of A';
        title(newtitle);
        
    case false
        %% No diagonal, instead plot time slices
        
        assert(exist('slice_vals', 'var') && ~isempty(time_vals), ...
            'If diagonal is not plotted but matrix slices instead, provide slice_vals = [y1,y2].');
        CI_lowbound_trial   = shiftdim(CItwosided(1, 1 + time_vals(1)/2, :));
        CI_highbound_trial  = shiftdim(CItwosided(2, 1 + time_vals(1)/2, :));
        CI_lowbound_rest    = shiftdim(CItwosided(1, 27 + time_vals(2)/2, :));
        CI_highbound_rest   = shiftdim(CItwosided(2, 27 + time_vals(2)/2, :));
        
        % values
        time_vals(1) = 1+ time_vals(1)/2;
        time_vals(2) = 27 + time_vals(2)/2;
        
        % plot CI fills
        hold on;
        xaxis_fill = [xaxis_values, fliplr(xaxis_values)];
        yaxis_fill_trial = [CI_highbound_trial; flipud(CI_lowbound_trial)];
        fillplot(1) = fill(xaxis_fill, yaxis_fill_trial, purple, 'facecolor',  ...
            purple, 'edgecolor', purple, 'facealpha', 0.15);
        yaxis_fill_rest = [CI_highbound_rest; flipud(CI_lowbound_rest)];
        fillplot(2) = fill(xaxis_fill, yaxis_fill_rest, 1, ...
            'facecolor',  pinkish, 'edgecolor', pinkish, 'facealpha', 0.15);
        
        % Info
        t_str = sprintf('Trial %is', slice_vals(1));
        r_str = sprintf('Rest %is', slice_vals(2));
        l_hdl = legend([fillplot(1), fillplot(2)], {t_str, r_str}, ...
            'Location', 'southoutside', 'AutoUpdate', 'off');
        try l_hdl.ItemTokenSize = [10 10 10]; end
        
        % plot lines
        plot(xaxis_values, result.mean(time_vals(1),:,:), 'Color', purple, 'LineWidth', 1)
        plot(xaxis_values, result.mean(time_vals(2),:,:), 'Color', pinkish, 'LineWidth', 1)
        yline(result.chancelevel + offset  + offset, 'k-'); % adds chance level line
        xline(reststart + offset  + offset, 'k--'); % adds end of task phase
        
        % significant values slices
        sig_vals1 = H_vals(time_vals(1), :);
        sig_vals2 = H_vals(time_vals(2), :);
        sig_pts1 = round(max(CI_highbound_trial)+4) * ones(size(CI_highbound_trial))';
        sig_pts1(~sig_vals1)    = NaN;
        sig_first1              = find(~isnan(sig_pts1), 1, 'first');
        sig_last1               = find(~isnan(sig_pts1), 1, 'last');
        dispv(1, '    %s timecourse first and last significant values: %i, %i.', ...
            t_str, sig_first1, sig_last1);
        plot(sig_pts1, 'Color', purple, 'LineWidth', 1.5); % significance line
        
        sig_pts2 = round(max(CI_highbound_trial)+3) * ones(size(CI_highbound_trial))';
        sig_pts2(~sig_vals2)    = NaN;
        sig_first2              = find(~isnan(sig_pts2), 1, 'first');
        sig_last2               = find(~isnan(sig_pts2), 1, 'last');
        dispv(1, '    %s timecourse first and last significant values: %i, %i.', ...
            r_str, sig_first2, sig_last2);
        plot(sig_pts2, 'Color', pinkish, 'LineWidth', 1.5); % significance line
        hold off;
        
        newtitle = 'E. Example training set test accs.';
        title(newtitle);
        
    otherwise
        error('Unknown case');
end

% ticks, labels, title
font_sz = 8;
axis tight; axis square;
xticks(customticks + offset  + offset); xticklabels(customlabels);
set(gca, 'XTickLabelRotation', 60); set(gca,'FontSize',font_sz);
ylabel(y_labels, 'Interpreter', 'none', 'FontSize', font_sz);
xlabel(x_labels, 'Interpreter', 'none', 'FontSize', font_sz);

end

%% OLD CODE FOR CALLING THIS FUNCTION

%% Corrected
%
% % Call
% fighdl2 = plot_timecourse(result, targetfile, group_resname, customticks, ...
%     customlabels, 1, 1);
%
% % Save filled plot
% dispv(1, '    Saving timecourse plot with CI color fill to target directory...');
% targetfile2 = [targetfile '_bonf_tc'];
% save_plot(targetfile2, fighdl2);
%
%
% %% Training data slices (only corrected for now)
%
% task_rest_tp = [10 10];
% fighdl3 = plot_timecourse(result, targetfile, group_resname, customticks, ...
%     customlabels, 1, 0, task_rest_tp);
%
% % Save
% targetfile3 = [targetfile '_slices_bonf'];
% save_plot(targetfile3, fighdl3);

%% Uncorrected

% do_uncorrected = 0; % disabled for now as they are not reported
%
% if do_uncorrected
%
%     % Call
%     fighdl1 = plot_timecourse(result, targetfile, group_resname, customticks, ...
%         customlabels, 0, 1);
%
%     % Save filled plot
%     try
%         dispv(1, '    Saving timecourse plot with CI color fill to target directory...');
%         targetfile1 = [targetfile '_uncorr_tc']; disp(['    Writing ' targetfile1 '.mat'])
%         save([targetfile1 '.mat'], 'result');
%         save_plot(targetfile1, fighdl1);
%         result.files.mats{end+1} = [targetfile1 '.mat'];
%         result.files.figures{end+1} = [targetfile1 '.fig'];
%     catch e
%         warning('Saving plot failed, see error...')
%         disp(['    Error: ' e.message]);
%     end
% end