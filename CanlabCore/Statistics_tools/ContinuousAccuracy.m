function acc = ContinuousAccuracy(obj, pattern, predrange, unit)
% Calculate forced choice accuracy for unit increases in continuous
% predictions.  Requires units to be ranked ordered from 1:end.  Will work
% for Gianaros or Pain datasets.  Accuracies are not penalized for missing
% cases.  Probably best to run this on single subjects and then aggregate
% accuracies across subjects.
%
% :Usage:
% ::
%
%     acc = ContinuousAccuracy(obj, pattern, unit)
%
% :Inputs:
%
%   **obj:**
%        fmri_data() object with data stacked by
%        increasing levels of prediction.  Make sure 
%        obj.Y includes the training labels
%
%   **pattern:**
%        fmri_data() object with weight pattern
%
%   **predrange:**
%        specify the range of predictions (e.g., 1:5)
%
%   **unit:**
%        specify the unit increase in prediction (e.g., 1 or 2)
%
% :Outputs:
%
%   **acc:**
%        accuracy of prediction for specified units
%
% :Examples:
% ::
%
%    acc = ContinuousAccuracy(dat, pine, 1:5, 1)
%
% ..
%    Original version: Copyright Luke Chang 12/2013
% ..

pexp = apply_mask(obj, pattern, 'pattern_expression', 'ignore_missing'); %Calculate pattern expression


%Create pairwise matrix and find positive values for lower triangle
%-need to force it to be 5 and fill in missing with NaNs
for i = predrange
    for j = predrange
        if ~any(obj.Y==i) || ~any(obj.Y==j)
            a(i,j) = nan;
        else
            a(i,j) = pexp(obj.Y == i) - pexp(obj.Y == j);
        end
    end
end
lt = double(tril(a) > 0);
lt(isnan(a)) = nan;

%remove by a factor of unit
i = 1;
while i < max(predrange)    
    lt(i:i+unit-1,i) = 0;
    i = i + 1;
end

n = length(predrange)-unit;
miss = tril(isnan(lt)); %calculate lower triangle accuracy
miss(logical(eye(length(predrange))))=0; %remove any nans on diagonal
total = (n*(n+1)/2) - sum(sum(miss));  %subtract out any nans from total
acc = sum(nansum(lt))/total;


