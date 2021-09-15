% [all_CONFMAT, confmat_fighdl] = afx_subject_timeresolved_confmats(cfg, results)
%
% Compiles confusion matrix data, creates a matrix with or without
%  undecided class (see TDT demo3 of TDT 3.999D or higher for an
%  explanation) and plots this data for a subject. Does so for all input
%  masks. 
%
% Can support an undecided class in case the classifier used gives
%  it as output in TDT (likely from TDT 4 onwards or an earlier beta 
%  version).
%
% Needs: results.confusion_matrix(_plus_undecided) from 
%   afx_tdt_subject_timeresolved_decoding.m
%  cfg with design files
%  optional: maskname
%
% Ingmar, 02-07-21

% Now supports a 4-th undecided class (classifier could not decide between 
% labels at all in majority vote) - Ingmar 12-07-21

function [all_CONFMAT, confmat_fighdl] = afx_subject_timeresolved_confmats(cfg, results)
all_CONFMAT = {}; CONFMAT = {}; confmat_fighdl = []; % init

%% Determine if an undecided class was added (compatible with TDT 3.999D)

confmat_str = 'confusion_matrix_plus_undecided';
if isfield(results, confmat_str) && any(strcmpi(cfg.results.output, confmat_str))
    dispv(1, '    Extracting confusion matrix values PLUS undecided values over time...');
    confmat_res = results.confusion_matrix_plus_undecided;
    undecided = 1;
elseif isfield(results, 'confusion_matrix') && any(strcmpi(cfg.results.output, 'confusion_matrix'))
    dispv(1, '    Extracting confusion matrix values over time...');
    confmat_res = results.confusion_matrix;
    undecided = 0;
else
    dispv(1, '  No confusion matrix was detected in result output, continuing to data saving...');
    return
end

%% Construct CONFMAT

% Enable this only to check if set IDs are ordered across
% time-resolved matrix correctly (note: remove the [0, 100] clim for the
% imagesc plotting below in plot_subject_confmats).
do_id_check = 0; % <--- keep 0 for normal analysis!

% Loop over masks
for m_ind = 1:length(cfg.files.mask)
    conf_mat_sz = size(confmat_res.set(1).output{m_ind});
    for row_ind = 1:conf_mat_sz(1)
        for col_ind = 1:conf_mat_sz(2)
            for set_ind = 1:length([confmat_res.set.set_id])
                
                % Get train x test timepoints
                curr_set    = confmat_res.set(set_ind).output{m_ind}; % get set
                curr_set_id = confmat_res.set(set_ind).set_id; % get set id
                curr_tp     = num2cell(cfg.design.set2conditions_map(curr_set_id, :));
                [~, train_ind, test_ind] = curr_tp{:}; 
                CONFMAT{row_ind, col_ind}.MAT(train_ind, test_ind) = curr_set(row_ind, col_ind);
                
                if do_id_check % Check order
                    CONFMAT_set_id{row_ind, col_ind}.MAT(train_ind, test_ind) = curr_set_id; %#ok<*UNRCH>
                end
            end
            
            % Check for nans (not sure if still needed)
            if any(isnan([CONFMAT{row_ind, col_ind}.MAT]))
                warning('Nans in CONFMAT{%i, %i}.output: %i nans in %i values', ...
                    row_ind, col_ind, sum(isnan([CONFMAT{row_ind, col_ind}.MAT])), ...
                    numel([CONFMAT{row_ind, col_ind}.MAT]));
            end
            
            if col_ind < 4
                % Get and set current labels
                label_val2str = table2cell(cfg.design.label_val2str_map);
                rowstr = label_val2str{row_ind, 2};
                colstr = label_val2str{col_ind, 2};
            elseif col_ind == 4
                rowstr = label_val2str{row_ind, 2};
                colstr = 'Undecided';
            end
            CONFMAT{row_ind, col_ind}.compstr = sprintf('True "%s" vs. predicted "%s"', rowstr, colstr); %#ok<*AGROW>
            CONFMAT{row_ind, col_ind}.elstr = sprintf('CONFMAT element (%0i, %0i)', row_ind, col_ind);
            dispv(1, '    Train-by-test time-points compiled for %s: %s', ...
                CONFMAT{row_ind, col_ind}.elstr, CONFMAT{row_ind, col_ind}.compstr);
            
            if do_id_check % Check order
                CONFMAT_set_id{row_ind, col_ind}.compstr = sprintf('set_id true "%s" vs. predicted "%s"', rowstr, colstr); %#ok<*AGROW>
                CONFMAT_set_id{row_ind, col_ind}.elstr = sprintf('CONFMAT set_id element (%0i, %0i)', row_ind, col_ind);
            end
        end
    end
    
    try
        [~, curr_maskname, ~] = fileparts(cfg.files.mask{m_ind});
    catch
        dispv(1,'    Getting maskname from file string failed.');
        curr_maskname = sprintf('mask-%02i', m_ind);
    end
        
    % Plotting confusion matrix timeseries for all elements x train x test
    % ------------------------------------------------------------ %
    confmat_fighdl = plot_confmats(cfg, CONFMAT, undecided, curr_maskname);
    % -------------------------------------------------------------%
    
    % Save plot -- filename; maskname has subject-id (at least for sim)
    plot_fname = fullfile(cfg.level2_dir, sprintf('sub-%02i', cfg.subj), ...
        sprintf('confusion-matrix_%s', curr_maskname));
    save_plot(plot_fname, confmat_fighdl, 1);
    
    if do_id_check  % Set_id plot
        setid_fighdl = plot_confmats(cfg, CONFMAT_set_id, undecided);
        plot_fname = fullfile(cfg.level2_dir, sprintf('sub-%02i', cfg.subj), ...
            sprintf('confusion-matrix_set-id_%s', curr_maskname));
        save_plot(plot_fname, setid_fighdl, 1);
    end
    
    % Store CONFMAT
    all_CONFMAT{m_ind} = CONFMAT;
end % masks

end % func




