% function [results, fighdl] = afx_plot_2d_mats_mean_ci95(data, labels, title_str, fighdl, sort_ind)
%
% short description
%
% IN
%  labels: either: cellstr, will use the same labels for x and y axis. or
%     struct with fields labels.x and labels.y, both with cellstr that will
%     be used for x and y respectively
% OPTIONAL
%  title_str, fighdl, sort_ind
% OUT
%  reults: struct with fields
%   .alphaval: value of alpha
%  fighdl to function to be used
%
% add as option: sort_ind
%
% TODO: Finish this description
%
% Kai 2020-11-11

function [results, fighdl] = afx_plot_2d_mats_mean_ci95(data, labels, title_str, fighdl, sort_ind)

%% defaults

if ~exist('title_str', 'var') || isempty(title_str)
    title_str = 'mean mat and ci95';
end

if ~exist('fighdl', 'var') || isempty(fighdl)
    fighdl = figure('name', title_str);
end

if ~exist('labels', 'var')
    labels ='';
end

if isempty(labels)
    labels.x = 1:size(data, 2);
    labels.y = 1:size(data, 1);
end


if ~isfield(labels, 'x')
    if isfield(labels, 'y')
        error('please provide both labels.x and .y')
    end
    org_labels = labels;
    clear labels
    labels.x = org_labels;
    labels.y = labels.x;
elseif ~isfield(labels, 'y')
    error('please provide both labels.x and .y')
end

assert(length(labels.x)==size(data, 2));
assert(length(labels.y)==size(data, 1));
%% resort data and labels if provided

if exist('sort_ind', 'var') && ~isempty(sort_ind)
    if length(unique(sort_ind)) ~= length(sort_ind)
        error('sort_ind contained at least one index multiple times. I guess that is not what you wanted.')
    end
    disp('Resorting matrix')
    if ~isfield(sort_ind, 'x')
        if isfield(sort_ind, 'y')
            error('please provide in sort_ind both .x and .y')
        end
        org_sort_ind = sort_ind;
        clear sort_ind
        sort_ind.x = org_sort_ind;
        sort_ind.y = sort_ind.x;
    else
        if ~isfield(sort_ind, 'y')
            error('please provide in sort_ind both .x and .y')
        end
    end
    
    data = data(sort_ind.y, sort_ind.x, :);
    labels.x = labels.x(sort_ind.x);
    labels.y = labels.y(sort_ind.y);
end

%% Calculate values
disp('Calculating values')
results.all_means = mean(data, 3);

results.all_counts = size(data, 3);

% results.all_stds = nanstd( data, 0, 3); % uses N-1 for std
% results.all_sems = nanstd( data, 0, 3) ./ sqrt(results.all_counts); % uses N-1 for std
results.all_CI95 = nanstd( 1.96 * data, 0, 3) ./ sqrt(results.all_counts); % uses N-1 for std

results.alphaval = 0.05;
[results.all_H, results.all_P, results.all_CI, results.all_STATS] = ttest(data, 0, 'dim', 3, 'alpha', results.alphaval);

%% Plot values

% mean
subplot(2, 3, 1)
imagesc(results.all_means)
mark_nans(results.all_means)

title(get(gcf, 'name'), 'interpreter', 'none');
colorbar;

xticks = 1:size(results.all_means, 2);
yticks = 1:size(results.all_means, 1);

set(gca, 'Xtick', xticks, 'XTickLabel', labels.x(xticks), 'XTickLabelRotation', 60);
set(gca, 'YTick', yticks, 'YTickLabel', labels.y(yticks));
axis equal
axis tight

% CI 95
subplot(2, 3, 2)
imagesc(results.all_CI95)
mark_nans(results.all_CI95)

title('CI95');
colorbar;

xticks = 1:size(results.all_means, 2);

set(gca, 'Xtick', xticks, 'XTickLabel', labels.x(xticks), 'XTickLabelRotation', 60);
axis equal
axis tight
% yticks = 1:size(results.all_means, 1);
% set(gca, 'YTick', yticks, 'YTickLabel', labels.y(yticks));


% ttest & p-values
subplot(2, 3, 3)
imagesc(results.all_P, [1e-7, 1])
set(gca, 'colorscale', 'log')
H_nan = results.all_H;
H_nan(results.all_H==0) = nan;
mark_nans(H_nan, 'x', 'k', '', 'k')

title(['ttest alpha=' num2str(results.alphaval)])
colorbar;
axis equal
axis tight

% mean vs CI95 (+ mean)
subplot(2, 3, 4)
plot(results.all_means(:), [results.all_means(:)-results.all_CI95(:), results.all_means(:)+results.all_CI95(:)], '.')
xlabel('mean')
ylabel('mean+-CI95')
% axis equal
axis tight
grid on

subplot(2, 3, 5)
plot(results.all_means(:), results.all_CI95(:), '.');
xlabel('mean(RHOs)')
ylabel('CI95(RHOs) abs value');
% axis equal
axis tight
grid on
