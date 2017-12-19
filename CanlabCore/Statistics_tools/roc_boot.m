function [ci, names] = roc_boot(input_vals, binary_outcome, thr, verbose)
% :Usage:
% ::
%
%     [ci, names] = roc_boot(input_vals, binary_outcome, thr, [verbose flag])
%
% Returns bootstrapped 95% confidence intervals for sensitivity,
% specificity, and PPV at a given threshold.
%
% :Examples:
% ::
%
%    thr = ROC.class_threshold
%    input_vals = input(ind);
%    binary_outcome = outcome(ind);


if nargin < 4, verbose = 1; end

bootfun = @(input_vals, binary_outcome) classification_stats(input_vals, binary_outcome, thr);

sens_spec_ppv = bootstrp(1000, bootfun, input_vals, binary_outcome);

% 95% CIs
for i = 1:3
    ci{i} = [prctile(sens_spec_ppv(:, i), 2.5)  prctile(sens_spec_ppv(:, i), 97.5)];
end

names = {'95% CI for sensitivity' '95% CI for specificity' '95% CI for PPV'};

if verbose
    ssp = bootfun(input_vals, binary_outcome);
    
    for i = 1:3
        fprintf('%s\t%.0f%% +- (%.0f-%.0f%%)\n', names{i}, ssp(i)*100, ci{i}(1)*100, ci{i}(2)*100);
    end
    
end


end % main function


function sens_spec_ppv = classification_stats(input_vals, binary_outcome, thr)

wh = input_vals >= thr;

tpr = sum(wh(binary_outcome)) ./ sum(binary_outcome);
fpr = sum(wh(~binary_outcome)) ./ sum(~binary_outcome);

npos = sum(tpr .* binary_outcome);
nfp = sum(fpr .* ~binary_outcome);
ppv = npos ./ (npos + nfp);

sens_spec_ppv = [tpr 1-fpr ppv];

end

