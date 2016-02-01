function [xvals, tpr, fpr, auc, c_bias] = roc_calc(input_vals, binary_outcome, xvals)
% Calculate Receiver Operating Characteristic plot (ROC) given P-values
%
% :Usage:
% ::
%
%     [xvals, tpr, fpr, auc, c_bias] = roc_calc(input_vals or input values, binary_outcome, [xvals : threshold vals to assess])
%
% :Inputs:
%
%   **input_vals:**
%        continuous-valued observations to classify (e.g., fMRI activity)
%
%   **binary_outcome:**
%        1 / 0 vector of which input observations are "hits"
%
%   **xvals:**
%        Criterion values you put in or every 10th percentile of the input
%        data distribution by default
%
% :Outputs:
%
%   **tpr:**
%        True positive rate for every step of ROC curve (sensitivity)
%
%   **fpr:**
%        False positive rate (1 - specificity)
%
%   **auc:**
%        Empirical estimate of area under the ROC curve
%
%   **c_bias:**
%        c measure of response bias at each step; MacMillan and Creelman 2005
%
% :Examples:
% ::
%
%    % May not work for p-values? may need to convert to t or something.
%    pvals = STATS.WTS.p;
%    isnull = DATA.true_weights == 0;
%    [xvals, tpr, fpr] = roc_calc(pvals, isnull);
%    figure; plot(fpr, tpr, 'ko-','Color', 'k', 'LineWidth', 2);
%
%    figure; plot(xvals, fpr, 'bo-')
%    hold on; plot([0 1], [0 1], 'k', 'LineWidth', 2);
%    set(gca, 'XLim', [0 .2], 'YLim', [0 .2])
%    xlabel('Nominal false positive rate');
%    ylabel('Actual false positive rate');
%
% :See Also: roc_plot.m

if ~islogical(binary_outcome), disp('Warning!! binary_outcome must be logical.'); end
binary_outcome = logical(binary_outcome);

if nargin < 3 || isempty(xvals)
    xvals = prctile(input_vals, [0:10:100]);
end

[tpr, fpr] = deal(zeros(size(xvals)));

indx = 1;
for x = xvals
    wh = input_vals >= x;
    
    tpr(indx) = sum(wh(binary_outcome)) ./ sum(binary_outcome);
    fpr(indx) = sum(wh(~binary_outcome)) ./ sum(~binary_outcome);
    
    indx = indx + 1;
end

auc = calc_auc(fpr, tpr);

c_bias = ( norminv(max(.0001, min(0.9999, tpr))) + norminv(max(.0001, min(0.9999, fpr))) ) ./ 2;

end % function





function auc = calc_auc(fpr, tpr)

[u, wh] = unique(fpr);
u2 = tpr(wh);

% fix for AUC = 1 if no overlap; triangle method not perfectly accurate
% here.
if any(u == 0 & u2 == 1), auc = 1; return, end

for i = 2:length(u)
    
    xdiff = u(i) - u(i - 1);
    ydiff = u2(i) - u2(i - 1);
    a(i) = xdiff * u2(i - 1) + xdiff * ydiff / 2;  % area of rect + area of triangle
    
end


auc = sum(a);

end





