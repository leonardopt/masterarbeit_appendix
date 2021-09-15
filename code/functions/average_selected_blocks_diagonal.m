% Function to average the diagonal of the selected blocks.

function [avmat, lowerCI, upperCI] = average_selected_blocks_diagonal(allmat, cond_n)

if cond_n > 3
    error('Function not implemented yet for more than three conditions')
end


% If congruent condition
if ~isnan(allmat{1}(1,1))
    
    alldiags = [];
    for m = 1:length(allmat)
        
        mat = allmat{m};
        aehdiag = diag(mat(1:45, 1:45));
        kondiag = diag(mat(46:90, 46:90));
        sumdiag = diag(mat(91:135,91:135));
        
        alldiags = vertcat(alldiags, [aehdiag kondiag sumdiag]');
        
    end 
    
    zmean = nanmean(atanh(alldiags), 1);
    zCI95 = 1.96 * ( nanstd( atanh(alldiags), 0, 1) ./ sqrt(size(alldiags, 1)));
    
    lowerCI = tanh(zmean-zCI95);
    upperCI = tanh(zmean+zCI95);
    
    avmat = tanh(zmean);
    
% If incongruent condition    
else
    
    alldiags = [];
    
    for m = 1:length(allmat)
        
        mat = allmat{m};
        diag1 = diag(mat(46:90, 1:45));
        diag2 = diag(mat(91:135, 1:45));
        diag3 = diag(mat(91:135, 46:90));
        diag4 = diag(mat(1:45, 46:90));
        diag5 = diag(mat(1:45, 91:135));
        diag6 = diag(mat(46:90, 91:135));
                
        alldiags = vertcat(alldiags, [diag1, diag2, diag3, diag4, diag5, diag6]');
        
    end 
    
    zmean = nanmean(atanh(alldiags), 1);
    zCI95 = 1.96 * ( nanstd( atanh(alldiags), 0, 1) ./ sqrt(size(alldiags, 1)));
    
    lowerCI = tanh(zmean-zCI95);
    upperCI = tanh(zmean+zCI95);
    
    avmat = tanh(zmean);
    
end

end