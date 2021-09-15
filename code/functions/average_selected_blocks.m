% Function to average the selected blocks.
% First, it divides the matrix in three parts of equal columns, which are
% then vertcat to one another  Result: matrix with dimensions 135*3 x 45

function [av_curve, lowerCI, upperCI] = average_selected_blocks(mat, cond_n)

% OLD FUNCTION
% if cond_n > 3
%     error('Function not implemented yet for more than three conditions')
% end 
% matcat = vertcat(mat(:,1:end/cond_n), mat(:,end/cond_n + 1:(end/cond_n) *2), mat(:,(end/cond_n) *2 +1:end));
% matfinal = matcat(all(~isnan(matcat),2),:);
% avmat = tanh(nanmean(atanh(matfinal)));
% 
% zCI95 = 1.96 * ( nanstd( atanh(matfinal), 0, 1) ./ sqrt(size(matfinal, 1)));
% 
% 
% zmean = nanmean(atanh(matfinal));
% 
% lowerCI = tanh(zmean-zCI95);
% upperCI = tanh(zmean+zCI95);

if cond_n > 3
    error('Function not implemented yet for more than three conditions')
end 

for i =1:length(mat)
    
    curr_mat = mat{i};

    matcat = vertcat(curr_mat(:,1:end/cond_n), curr_mat(:,end/cond_n + 1:(end/cond_n) *2), curr_mat(:,(end/cond_n) *2 +1:end));
    matfinal = matcat(all(~isnan(matcat),2),:);
    avmat = tanh(nanmean(atanh(matfinal)));
    zallmeans(i, :) = nanmean(atanh(matfinal));

end




% Mean over all participants 

zmean = nanmean(zallmeans, 1);

zCI95 = 1.96 * ( nanstd( zallmeans, 0, 1) ./ sqrt(size(zallmeans, 1)));



lowerCI = tanh(zmean-zCI95);
upperCI = tanh(zmean+zCI95);

av_curve = tanh(zmean);


end

