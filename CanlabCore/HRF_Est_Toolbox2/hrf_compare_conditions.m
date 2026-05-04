function cmp = hrf_compare_conditions(condA, condB, varargin)
%HRF_COMPARE_CONDITIONS Compare two condition-averaged trial sets.
% cmp = hrf_compare_conditions(condA, condB)
% condA/condB are outputs from hrf_average_condition_trials

p = inputParser;
p.addRequired('condA', @isstruct);
p.addRequired('condB', @isstruct);
p.addParameter('Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.parse(condA, condB, varargin{:});
opts = p.Results;

if numel(condA.time) ~= numel(condB.time)
    error('Condition windows must have same number of time points.');
end

meanA = condA.mean(:);
meanB = condB.mean(:);
diffMean = meanA - meanB;

% Summary features
[peakA, iA] = max(meanA);
[peakB, iB] = max(meanB);
aucA = trapz(condA.time, meanA);
aucB = trapz(condB.time, meanB);

% Timepoint-wise statistics (if stats toolbox exists)
pvals = nan(size(meanA));
tvals = nan(size(meanA));
for t = 1:numel(meanA)
    try
        [~, p, ~, stats] = ttest2(condA.trials(:, t), condB.trials(:, t), 'Vartype', 'unequal');
        pvals(t) = p;
        tvals(t) = stats.tstat;
    catch
        % Leave NaN if ttest2 is unavailable
    end
end

cmp = struct();
cmp.time = condA.time(:);
cmp.conditionA = condA.condition;
cmp.conditionB = condB.condition;
cmp.meanA = meanA;
cmp.meanB = meanB;
cmp.mean_diff = diffMean;
cmp.peak = struct('A', peakA, 'A_time', condA.time(iA), 'B', peakB, 'B_time', condB.time(iB));
cmp.auc = struct('A', aucA, 'B', aucB, 'diff', aucA - aucB);
cmp.p_value = pvals;
cmp.t_value = tvals;
cmp.significant = pvals < opts.Alpha;
cmp.alpha = opts.Alpha;
end
