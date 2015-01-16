function statsimg = ttest(fmridat, pvalthreshold, thresh_type)
% T-test on fmri_data class object
% statsimg = ttest(fmridat, pvalthreshold, thresh_type)
%
% ttest(fmridat, p-value threshold, thresh_type)
%
% p-value threshold: p-value, e.g., .05 or .001 or [.001 .01 .05]
% thresh_type: 'uncorrected', 'fwe', or 'fdr'
%
% e.g., 
% T-test, Construct a stats_image object, threshold and display:
% statsimg = ttest(fmridat, .001, 'unc');
% orthviews(statsimg);
%
% Re-threshold and display:
% statsimg = threshold(statsimg, .000001, 'unc');
% orthviews(statsimg);
% 
% statsimg = threshold(statsimg, .01, 'fdr');
% orthviews(statsimg);
%
% NOTE: for two-sample T-test, use fmri_data.regress

fprintf('One-sample t-test\n')
fprintf('Calculating t-statistics and p-values\n');

err = ste(fmridat.dat')';

t = nanmean(fmridat.dat')' ./ err;

N = single(sum(~(isnan(fmridat.dat') | fmridat.dat' == 0) , 1));

df = size(fmridat.dat, 2) - 1;

% two-tailed - now done in statistic_image constructor
% p = 2 * (1 - tcdf(abs(t), df));

% t = t';
% p = p';

% removed: done in constructor: 'p', p, 'p_type', 'two-tailed',
statsimg = statistic_image('type', 't', 'dat', t, 'ste', err, 'N', N, 'volInfo', fmridat.mask.volInfo, 'dfe', df);

% get rid of NaNs for empty vals
stats.p(isnan(t)) = 1;
stats.t(isnan(t)) = 0;

statsimg.removed_voxels = fmridat.removed_voxels;

if nargin > 2
    
    statsimg = threshold(statsimg, pvalthreshold, thresh_type);
    
end



end
