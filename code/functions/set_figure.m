% Function to set figure name depending on data type (whole brain,
% simulated, ROI)


function set_figure(cfg, name_fig, position)

data_type  = cfg.data_type;


if isfield(cfg, 'voxelrange')
    voxelrange = cfg.voxelrange;
    mean_n_voxels = cfg.mean_n_voxels;
    if numel(voxelrange) == 1
        voxelrange = num2str(voxelrange);
    else
        minv = num2str(voxelrange(1));
        maxv = num2str(voxelrange(2));
        voxelinfo = sprintf('n_voxels[min=%i,max=%i,mean=%g]', voxelrange(1), voxelrange(2), round(mean_n_voxels));
    end
end

switch data_type
    
    case 'wholebrain'
        if exist('voxelrange', 'var')
            figure('name', sprintf('%s', name_fig), ...
                'Position',  position);
            title(sprintf('%s\nwhole brain -  %s', name_fig, voxelinfo))
        else
            figure('name', sprintf('%s', name_fig), ...
                'Position',  position);
            title(sprintf('%s\nwhole brain', name_fig))
        end
        
    case 'simulated'
        source_data = cfg.sim_name;
        figure('name', sprintf('%s', name_fig), ...
            'Position',  position);
        title(sprintf('%s\n%s: [%s] voxels', name_fig, source_data, voxelrange(1)))
        
    case 'ROI'
        if exist('voxelrange', 'var')
            
            roiname    = cfg.ROI;
            figure('name', sprintf('%s', name_fig), ...
                'Position',  position);
            title(sprintf('%s\nROI(%s)-  %s', name_fig, roiname, voxelinfo))
        else
            roiname    = cfg.ROI;
            figure('name', sprintf('%s', name_fig), ...
                'Position',  position);
            title(sprintf('%s\nROI_%s', name_fig, roiname))
            
        end
        
    otherwise
        error('Error in setting figure. Please check input')
        
        
        
end