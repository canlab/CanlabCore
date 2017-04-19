function [oddeven_r_spearmanbrown, tabledat] = reliability(DAT, varname, varargin)

% DAT = study_canlab_dataset{i};
% varname = 'NPS_cosine_similarity_none';
% optional_inputs = { 'conditional', {'ok_trials' 1} };

% spearman-brown prophecy formula
spearmanbrown = @(r) 2*r / (1+(2-1)*r);

[~, celldat] = get_var(DAT, varname, varargin{:});

clear means robmeans means12

for j = 1:length(celldat)
    % For each subject
    
    t = length(celldat{j});
    
    % odd/even reliability
    wh1 = 1:2:t;
    wh2 = 2:2:t;
    
    means(j, 1) = nanmean(celldat{j}(wh1));
    means(j, 2) = nanmean(celldat{j}(wh2));
    
    robmeans(j, 1) = robust_mean(celldat{j}(wh1));
    robmeans(j, 2) = robust_mean(celldat{j}(wh2));
    
    % first/second half reliability
    wh1 = 1:floor(t/2);
    wh2 = floor(t/2):t;
    
    means12(j, 1) = nanmean(celldat{j}(wh1));
    means12(j, 2) = nanmean(celldat{j}(wh2));
    
    
end

mean_trials = mean(cellfun(@length, celldat));

[oddeven_r, oddeven_p] = corr(means(:, 1), means(:, 2));
[oddeven_robustmean_r, oddeven_robustmean_p] = corr(robmeans(:, 1), robmeans(:, 2));

oddeven_r_spearmanbrown = spearmanbrown(oddeven_r);
oddeven_r_robustmean_spearmanbrown = spearmanbrown(oddeven_robustmean_r);

[firstsecondhalf_r, firstsecondhalf_p] = corr(means12(:, 1), means12(:, 2));
%[firstsecondhalf_robustmean_r, firstsecondhalf_robustmean_p] = corr(robmeans(:, 1), robmeans(:, 2));

firstsecondhalf_r_spearmanbrown = spearmanbrown(firstsecondhalf_r);
%firstsecondhalf_r_robustmean_spearmanbrown = spearmanbrown(firstsecondhalf_robustmean_r);

tabledat = table(mean_trials, oddeven_r_spearmanbrown, oddeven_p, oddeven_r_robustmean_spearmanbrown, ...
    oddeven_robustmean_p, firstsecondhalf_r_spearmanbrown, firstsecondhalf_p);

end % function
