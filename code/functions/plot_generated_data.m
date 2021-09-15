% cfg = plot_generated_data(cfg, subj, ts_data, weighted_data, nii_dat, nii_nam, offset)
%
% Plots several aspects of the generated data (timeseries, voxels, etc.)
% and saves them to pre-specified locations
%
% IN: ts_data (timeseries data), nii_dat (nifti run-data after writing),
%       nii_nam (nifti names), cfg.saa.data_gen
%
% NEEDS: save_plot.m (which can use TDT)
%
% Ingmar, 14-08-20

% Added timecourse heatmap, Ingmar 08-01-21

function cfg = plot_generated_data(cfg, subj, ts_data, weighted_data, ...
    all_nii_dat, nii_nam, offset)
%% Turn plots on/off

do.timeseries       = 1;
do.timeseries_noise = 1;
do.single_voxel_ts  = 1;

%% Get parameters

layer_number    = cfg.saa.data_gen.layer_number;
layer_descr     = cfg.saa.data_gen.layer_descr;
layertotal      = length(layer_number);
aftereffect     = cfg.saa.data_gen.aftereffect;

%% Checks

% If plotting is needed
if isfield(cfg.saa.data_gen, 'plotting') && ~cfg.saa.data_gen.plotting
    dispv(1, '    Skipping subject plots since cfg.saa.data_gen.plotting == 0...\n');
    return % <--- Skip function
end

% If vectorized graphics will be saved
if ~isfield(cfg.saa.data_gen, 'save_svg')
    cfg.saa.data_gen.save_svg = 1; % auto-save
end

%% Target directory

targetdir_sub = fullfile(cfg.preproc_dir, sprintf('sub-%02i', subj));

%% Create timeseries fig

if do.timeseries
    dispv(1, '  A1. Visually check timeseries after model generation (first run)...');
    
    tsplotname = sprintf('Raw data timeseries (post-convolution) for sub-%02i', subj);
    tsfig      = figure('name', tsplotname, 'Position', [800 400 800 200*layertotal]);
    
    for layer_ind = 1:layertotal
        % Create subplot
        if layertotal <= 4 || mod(layertotal, 2) ~= 0 % uneven number
            subplot(layertotal, 1, layer_ind);
            ylabel(string(layer_ind), 'Interpreter', 'none');
            set(gcf, 'Position', [800 400 800 200*layertotal]);
        else
            subplot(layertotal/2, 2, layer_ind);
            set(gcf, 'Position', [400 500 900 150*(layertotal/2)]);
            title(sprintf('Set %i', layer_ind), 'Interpreter', 'none'); % add snr to subtitle
        end
        
        % Add model name if available (e.g. for plot name)
        try
            curr_aftereffect = replace(aftereffect, 'sustainedsuppression', 'suppr.');
            curr_modelset_name = sprintf('Set %i (%s, %s)', layer_ind, curr_aftereffect, ...
                layer_descr{layer_ind});
            cfg.saa.data_gen.fullmodelname = curr_modelset_name; % store
            dispv(2, '    Creating plot for: %s...', curr_modelset_name);
        catch
            dispv(1, ['    Could not create name for current timeseries model: ', ...
                'cfg.saa.data_gen.layer_descr is likely missing']);
            cfg.saa.data_gen.fullmodelname = string(layer_ind); % store
        end
        
        % Create subplot
        data_to_plot = ts_data{layer_ind}.X; % <--- Change dimension here for different run/block
        plot(data_to_plot);
        
        title(cfg.saa.data_gen.fullmodelname, 'Interpreter', 'none');
        axis tight;  xticks(round(linspace(0, size(data_to_plot, 1), 9)));
        curr_ylim = get(gca, 'YLim'); ylim([curr_ylim(1)*1.3 curr_ylim(2)*1.1]); % extend fig bounds
    end % layers
    
    suptitle('Timeseries (TR) convolved with HRF');
    
    % Save fig
    ts_targetfile = fullfile(targetdir_sub, ['timeseries-plot_' tsplotname]); % now adds plot title
    ts_targetfile = lower(strrep(strrep(strrep(ts_targetfile, ' ', '_'), ...
        '(', ''), ')', '')); % very ugly but seems there's no better way?
    
    dispv(1, '    Saving timeseries plot "%s" to simulation directory...', tsplotname);
    if ~cfg.saa.data_gen.save_svg
        dispv(1,'   "%s" was not saved as .svg because cfg.saa.data_gen.save_svg == 0', tsplotname);
        saveas(tsfig, [ts_targetfile '.png']); % Add .png otherwise saved as .fig
    else
        save_plot(ts_targetfile, tsfig, 1) % no fig
    end
    
    %% NEW: Add heatmap plot of all voxels side by side
    
    dispv(1, '  A2. Visually check weighted data...');
    
    tsplotname2 = sprintf('Weighted data timeseries per parameter set, 7th run, sub %02i', subj);
    tsfig2      = figure('name', tsplotname2, 'Position', [400 800 800 250*layertotal]);
    
    for layer_ind = 1:layertotal
        % Create subplot
        if layertotal <= 4 || mod(layertotal, 2) ~= 0 % uneven number
            subplot(layertotal, 1, layer_ind);
            ylabel(string(layer_ind), 'Interpreter', 'none');
        else
            subplot(layertotal/2, 2, layer_ind);
            set(gcf, 'Position', [400 800 800 200*(layertotal/2)]);
            title(sprintf('Set %i', layer_ind), 'Interpreter', 'none'); % add snr to subtitle
        end
        
        data_to_plot = weighted_data{layer_ind}(7).X' - offset; % zero-center using SPM offset
        imagesc(data_to_plot); colormap(customcolormap_preset('pasteljet'));
        axis tight; colorbar; caxis([min(min(data_to_plot)) max(max(data_to_plot))])
    end % layers
    
    suptitle(tsplotname2);
    
    % Save fig
    ts_targetfile2 = fullfile(targetdir_sub, ['plot_' tsplotname2]); % now adds plot title
    ts_targetfile2 = lower(strrep(strrep(strrep(ts_targetfile2, ' ', '_'), ...
        '(', ''), ')', '')); % very ugly but seems there's no better way?
    
    dispv(1, '    Saving timeseries plot "%s" to simulation directory...', tsplotname2);
    if ~cfg.saa.data_gen.save_svg
        dispv(1,'   "%s" was not saved as .svg because cfg.saa.data_gen.save_svg == 0', tsplotname2);
        saveas(tsfig2, [ts_targetfile2 '.png']); % Add .png otherwise saved as .fig
    else
        save_plot(ts_targetfile2, tsfig2, 1) % no fig
    end
    
    %%  NEW: TS and stick function combined
    
    if contains(cfg.saa.data_gen.aftereffect, 'suppression') || contains(cfg.saa.data_gen.aftereffect, 'fadeoff')
        dispv(1, '  A3. Signal model for first block stick function and TS (for showing in paper)');
        ts3plotname = sprintf(' sub-%02i', subj);
        ts3fig      = figure('name', ts3plotname, 'Position', [400 400 400 400]);
        
        hold on;
        data = weighted_data{1}.X; data = data - offset;
        plot(data(1:60, :), 'Color', [237/255 111/255 164/255]); % plot first block
        
        % x-axis: stick onsets
        sticks  = [0 (1.25:3.25:27) + 1 (27:3.25:45)]; % task offset of 1 TR
        
        % lazy way (downsampling from models, which are in dt, does not work well)
        downscale = 3.5; y = cfg.saa.data_gen.response(2); % 8
        if contains(cfg.saa.data_gen.aftereffect, 'suppression')
            au = [1 y y y y y y y y 0 -4.7995 -2.7372 -1.4082 -0.5518 0] ./ downscale; % all stick AU
        elseif contains(cfg.saa.data_gen.aftereffect, 'fadeoff')
            au = [1 y y y y y y y y 0 -4.7995 -2.7372 -1.4082 -0.5518 0] ./ downscale; % all stick AU
        end
        for au_ind = 1:length(au)
            plot([sticks(au_ind)-0.001 sticks(au_ind) sticks(au_ind)+0.001], [0 au(au_ind) 0], 'k-');
        end
        yline(0, 'k-')
        
        hold off;
        
        curr_aftereffect = replace(aftereffect, {'sustainedsuppression', 'fadeoff', 'none'}, ...
            {'A. Active inhibition model', 'B. Fade-off model', 'C. No aftereffect'});
        title(curr_aftereffect, 'Interpreter', 'none');
        axis square; axis tight;
        xticks(0:10:60);  yticks(-2:0:2); ylabel('A.U.');
        curr_ylim = get(gca, 'YLim'); ylim([curr_ylim(1)*1.2 curr_ylim(2)*1.2])
        
        %  Save fig
        ts3_targetfile = fullfile(targetdir_sub, ['timeseries-plot_' ts3plotname]); % now adds plot title
        ts3_targetfile = lower(strrep(strrep(strrep(ts3_targetfile, ' ', '_'), ...
            '(', ''), ')', '')); % very ugly
        
        dispv(1, '    Saving timeseries plot "%s" to simulation directory...', ts3plotname);
        save_plot(ts3_targetfile, ts3fig, 1) % no fig
    end
    
end % do.timeseries

%% Create noise profile over time

if do.timeseries_noise
    dispv(1, '  B. Noise distribution over time for first block, first run...');
    
    noiseplotname   = sprintf('Noise profile (post-convolution) for sub-%02i', subj);
    noisetsfig      = figure('name', noiseplotname, 'Position', [800 400 600 150*layertotal]);
    
    % assume first layer has no (or least) noise
    timespan    = 1:60;
    refmodel    = weighted_data{1}(1).X(timespan, :) - offset;
    
    for layer_ind = 1:length(layer_number)
        
        % Get mean of all cues, tasks, and rests across run
        curr_cue_noise  = weighted_data{layer_ind}(1).X(timespan, 1:3) - offset;
        curr_task_noise = weighted_data{layer_ind}(1).X(timespan, 4:6) - offset;
        curr_rest_noise = weighted_data{layer_ind}(1).X(timespan, 7:9) - offset;
        
        % Create subplot
        if layertotal <= 4 || mod(layertotal, 2) ~= 0 % uneven number
            subplot(layertotal, 1, layer_ind);
        else
            subplot(layertotal/2, 2, layer_ind);
            set(gcf, 'Position', [800 100 833 815]);
        end
        
        hold on;
        plot(curr_cue_noise, 'c.', 'LineWidth', 0.5);
        plot(curr_task_noise, 'm.', 'LineWidth', 0.5);
        plot(curr_rest_noise, '.', 'Color', [1 0.8 0], 'LineWidth', 0.5); % more orange than yellow for visibility
        plot(refmodel, 'Color', [0.6 0.6 0.6]); % plot model over scatter
        hold off;
        
        % Info
        if layertotal <= 4 || mod(layertotal, 2) ~= 0 % uneven number
            subplot(layertotal, 1, layer_ind);
        else
            subplot(layertotal/2, 2, layer_ind);
            set(gcf, 'Position', [800 100 833 815]);
        end
        mean_SNR = round(cfg.saa.snr_db_total(layer_ind), 2);
        title(sprintf('Set %i: total SNR = %.2f dB', ...
            layer_ind, mean_SNR), 'Interpreter', 'none'); % add snr to subtitle
        ylabel('A.U.'); axis tight;  xticks(round(linspace(0, size(timespan, 2), 9)));
        curr_ylim = get(gca, 'YLim'); ylim([curr_ylim(1)*1.1 curr_ylim(2)*1.1]); % extend figure bounds
    end % layers
    
    currtitle = sprintf('Post-HRF noise scatter per param. set (1st run, %i block, %i TRs, offset removed)', ...
        round(length(timespan)/45), length(timespan));
    suptitle(currtitle);
    
    %% Save fig
    noise_targetfile = fullfile(targetdir_sub, ['timeseries-plot_' noiseplotname]); % now adds plot title
    noise_targetfile = lower(strrep(strrep(strrep(noise_targetfile, ' ', '_'), ...
        '(', ''), ')', '')); % very ugly
    
    dispv(1, '    Saving timeseries plot "%s" to simulation directory...', noiseplotname);
    if isfile(noise_targetfile), warning('Overwriting...'); end % warn
    if ~cfg.saa.data_gen.save_svg
        dispv(1,'   "%s" was not saved as .svg because cfg.saa.data_gen.save_svg == 0', noiseplotname);
        saveas(noisetsfig, [noise_targetfile '.png']); % Add .png otherwise saved as .fig
    else
        save_plot(noise_targetfile, noisetsfig, 1) % no fig
    end
end % do.timeseries_noise

%% Plot differences in voxel layer timeseries

if do.single_voxel_ts && ~isempty(all_nii_dat)
    dispv(1, '  C. Visually check differences in voxel layer timeseries after nifti writing...');
    
    mixed_count = 0; nii_dat = all_nii_dat{1};
    
    for row_ind = 1:size(nii_nam, 1)
        for col_ind = 1:size(nii_nam, 2)
            
            % Check n-th iteration of mixed effect voxels to skip (avoid
            % multiple mixed effect plots)
            try
                split_nii_nam = strsplit(string(nii_nam(row_ind, col_ind, 1, 1)), '_');
                if strcmpi(split_nii_nam(1), 'mixed')
                    mixed_count = mixed_count + 1;
                end
            catch e
                dispv(1, '    Naming voxel (%i, %i) failed: %s', row_ind, col_ind, e.message)
                continue;
            end
            if mixed_count > 2, break; end
            
            % Construct subplot
            diff_fig = figure('name', 'Single-voxel timeseries for first run', ...
                'Position', [800 400 400 700]);
            
            for p_ind = 1:size(nii_dat, 3)
                
                % Create subplot
                subplot(layertotal, 1, p_ind); % Plot all layers for 1 voxel
                plot(squeeze(nii_dat(row_ind, col_ind, p_ind, :)), 'k-') % after nifti creation
                curr_title = strrep(string(nii_nam(row_ind, col_ind, p_ind, 1)), '_', ' ');
                curr_title = strrep(curr_title, 'fullmodel', '_');
                title(sprintf('%s. %s', char(64 + p_ind), curr_title), ...
                    'Interpreter', 'none');
                axis tight;
            end % subplot
            
            %% Save plot
            
            vox_targetfile = fullfile(targetdir_sub, ['voxel-ts_' char(curr_title)]); % now adds plot title
            vox_targetfile = lower(strrep(vox_targetfile, ' ', '_'));
            
            if ~cfg.saa.data_gen.save_svg
                dispv(1,'   Plot was not saved as .svg because cfg.saa.data_gen.save_svg == 0');
                saveas(diff_fig, [vox_targetfile '.png']);
            else
                save_plot(vox_targetfile, diff_fig, 1)
            end
        end % col
    end % row
end % do

%% Done, close all plots

pause(3); close ALL;
