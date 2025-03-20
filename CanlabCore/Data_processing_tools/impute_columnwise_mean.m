% Impute Column-wise Mean
function [X, wh_good, wh_bad] = impute_columnwise_mean(X, nan_values)
% For each column of X, replace values matching any of the list in nan_vales with the non-nan-valued mean of the column
%
% X = impute_columnwise_mean(X, nan_values)
% 
% 

for i = 1:size(X, 2)

    for j = 1:length(nan_values)      % identify bad cases
        wh_bad(:, j) = X(:, i) == nan_values(j);
    end
    wh_good = ~any(wh_bad, 2);           % which cases have good values
    
    m(i) = mean(X(wh_good, i));          % mean of good cases

    X(~wh_good, i) = m(i);               % impute
end

end % main function
