function [se_within, stats] = barplot_get_within_ste(dat,varargin)
% :Usage:
% ::
%
%     [se_within, stats] = barplot_get_within_ste(dat)
%
% Compute within-subjects standard errors 
%
% :Note: The old version of this used average s.e.'s for contrasts
% of interest, which depend on contrast scaling and are less appropriate
% for use as error bars.
%
% The new version as of 12/10 uses the Loftus & Masson 1994 method and has
% been checked against the data in that paper.
%
% See help below for the L & M data and additional
% code for getting the ANOVA decomposition.
%
% Useful for barplots of conditions when calculating within-subject SEs
%
% :Examples:
% ::
%
%    %% data from Loftus & Masson, Table 2
%
%    mtx = [1 10 13 13 12.00
%    2 6 8 8 7.33
%    3 11 14 14 13.00
%    4 22 23 25 23.33
%    5 16 18 20 18.00
%    6 15 17 17 16.33
%    7 1 1 4 2.00
%    8 12 15 17 14.67
%    9 9 12 12 11.00
%    10 8 9 12 9.67];
%
%    dat = mtx(:, 2:4);  % data from Loftus & Masson, Table 2
%
%    [se_within, stats] = barplot_get_within_ste(dat)
%    fprintf('Within ste: %3.2f, 95%% CI: mean +/- %3.2f\n', se_within, stats.ci);
%
% Extra stuff from ANOVA table
%
% Mean square for condition: Variance of condition means * sample
% size...average squared variance accounted for by condition means
% ::
%
%    [n, k] = size(dat); 
%    MS_cond = var( mean(dat) - mean(dat(:)) ) * n;
%
%    MS_subject = var( mean(dat') - mean(dat(:)) ) * k;
%
%    datv = dat(:);
%    MS_total = scale(datv, 1)' * scale(datv, 1);

% notes: 6/3/2018 - tor added wh_omit line to omit all-NaN cases and
% row-wise replacement to make error bars robust to presence of NaNs in
% data.  Previously returned NaN/zero error bars in presence of NaNs

if nargin > 1 && ~isempty(varargin{1})
    
    covs = varargin{1};

    if size(covs,1) ~= size(dat,1), error('Covs and dat must have same no. of observations.'); end
    
    disp('Warning: covariates are no longer removed for within-subject SE bars.');
    
end

% remove cases with all NaN values
% otherwise, these cases could influence calculations and interfere with
% scaling below.
wh_omit = all(isnan(dat), 2); dat(wh_omit, :) = []; 

% Replace remaining NaNs with row means 
m = nanmean(dat')';
m = repmat(m, 1, size(dat, 2));
wh = isnan(dat);
dat(wh) = m(wh);

% number of subjects and conditions
[n, k] = size(dat); 

df_s = n - 1;  % degrees of freedom for subjects
df_c = k - 1;  % degrees of freedom for conditions
df_sxc = df_s * df_c;  % degrees of freedom for subj x condition

% within-subject, within-condition errors
errs = scale(scale(dat, 1)', 1)';

% Sums of squared errors
SS_sxc = errs(:)' * errs(:);

% Mean squared error, accounting for reduced df due to estimating marginal means
MS_sxc  = SS_sxc / df_sxc;

se_within = sqrt(MS_sxc / n);

ci = se_within * tinv(1 - (.05 / 2), df_sxc);

stats.se_within = se_within;
stats.ci = ci;
stats.ci_descrip = '95% confidence interval';
stats.SS_sxc = SS_sxc;
stats.MS_sxc = MS_sxc;

% Extra stuff from ANOVA table
% ---------------------------------------------------------------
% Mean square for condition: Variance of condition means * sample
% size...average squared variance accounted for by condition means
stats.MS_cond = var( mean(dat) - mean(dat(:)) ) * n;

stats.MS_subject = var( mean(dat') - mean(dat(:)) ) * k;

datv = dat(:);
stats.MS_total = scale(datv, 1)' * scale(datv, 1);

return



