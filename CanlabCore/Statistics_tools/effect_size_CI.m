function [d_CI, P, P_CI] = effect_size_CI(d, n, varargin)
% Confidence interval for effect size (Cohen's d or Hedges' g) using noncentral t distribution
%
% d_CI = effect_size_CI(d, n, [dfe])
%
% - This function assumes dfe is n-1 unless specified otherwise
% - Assumes normally distributed data
%
% Outputs:
% d_CI = 95% confidence interval for d, 2-tailed
% P = P-value for the input d
% P_CI = 95% confidence interval for P, 2-tailed
%
% :Inputs:
%
%   **d:**
%        Effect size estimate
%
%   **n:**
%        Sample size
%
% :Optional Inputs:
%   **'method', 'noncentralt':**
%        Use noncentral t method instead of Hedges and Olkin
%
%   **'n2', [integer]:**
%        Sample size for 2nd group
%        In this case, effect_size_CI returns the CI for the two-sample
%        effect size.  n is the first group's sample size.
%
%   **'alphaval', [numeric value]:**
%        Alpha value, acceptable false positive rate
%
%   **'tails', [1 or 2]:**
%        1 for one-tailed, 2 for two-tailed tests/CIs. Default = 2
%
% :Outputs:
%
%   **d_CI:**
%        95% confidence interval for d, 2-tailed by default
%
%   **P:**
%        P-value for the input d
%
%   **P_CI:**
%        P-values for lower and upper 95% confidence interval bounds on d, 2-tailed by default
%
%
% :Notes:
% ::
% --------------------------------------------------
%
% There are multiple methods for this and formulas are suggested by
% Hedges (1981), Theory for Glass's Estimator of Effect Size and Related Estimators
%
% This function implements two methods:
% 1.Hedges-Olkin, a parametric calculation that is quick to implement
% - the default
%
% Hedges LV, Olkin I. Statistical methods for meta-analysis. Orlando: Academic Press Inc; 2014. p. 86.
% https://www.statisticshowto.datasciencecentral.com/hedges-g/
% https://books.google.com/books?id=5obZnfK5pbsC&pg=PA10&dq=hedges+g&hl=en&sa=X&ved=0ahUKEwjLl6eL4NPRAhVowYMKHQzGDpsQ6AEIGjAA#v=onepage&q=hedges%20g&f=false
%
% 2. A version based on the non-central t-distribution
% - this is an "alpha version" and may not be correct
% the method suggested by David C. Howell, University of Vermont
% Using the CDF of the non-central t distribution for the observed t score
% https://www.uvm.edu/~statdhtx/methods8/Supplements/MISC/Confidence%20Intervals%20on%20Effects/Effect%20Sizes%20Confidence%20Intervals.html
%
% NOTE: Confidence intervals are not symmetric.
%
% :See also:
% ::
% --------------------------------------------------
% power_calc.m
%
% Formulas
% p2t = @(p, dfe, tails) abs(tinv(p ./ tails, dfe));
% p2d = @(p, dfmodel, dfe, tails) abs(tinv(p ./ tails, dfe)) ./ sqrt(dfmodel + dfe);  % dfmodel + dfe = N observations
% d2t = @(d, dfmodel, dfe) d .* sqrt(dfmodel+dfe);
% d2p = @(d, dfmodel, dfe, tails) tails .* (1 - tcdf(d .* sqrt(dfmodel+dfe), dfe));
% d2p = @(d, n, dfe, tails) tails .* (1 - tcdf(d .* sqrt(n), dfe));
%
%
% :Examples:
% ::
% ---------------------------------------
% e.g.,
% n = 100; d = 0.65;
% d_CI = effect_size_CI(d, n);
%
% [d_CI, P, P_CI] = effect_size_CI(0.5, 30)
%
% A two-sample version of the above (wider confidence interval):
% [d_CI, P, P_CI] = effect_size_CI(0.5, 30, 'n2', 30)
%
% % Generate 1 and 2-sample values for a range of effect sizes and plot
% % ---------------------------------------
% dvals = [0.01:0.05:2]';
% n = 20;
% [d_CI1, d_CI2] = deal(zeros(length(dvals), 2));
% P = [];
% 
% for i = 1:length(dvals)
% 
%     [d_CI1(i, :), P(i, 1)] = effect_size_CI(dvals(i), n);
%  
%     d_CI2(i, :) = effect_size_CI(dvals(i), n, 'n2', n);
% 
% end
% 
% figure; 
% subplot(1, 2, 1);
% plot(dvals, dvals, 'b', 'MarkerFaceColor', [.3 0 1]);
% hold on; 
% han = errorbar(dvals, dvals, dvals - d_CI1(:, 1), d_CI1(:, 2) - dvals);
% set(han, 'Color', 'b', 'LineWidth', 0.1)
% 
% han = errorbar(dvals(P < 0.05), dvals(P < 0.05), dvals(P < 0.05) - d_CI1(P < 0.05, 1), d_CI1(P < 0.05, 2) - dvals(P < 0.05));
% han2 = plot_horizontal_line(0); set(han2, 'Linestyle', '--');
% xlabel('Effect size (d)');
% ylabel('Effect size (d) with CIs')
% title('One-sample test');
% 
% subplot(1, 2, 2);
% plot(dvals, dvals, 'b', 'MarkerFaceColor', [.3 0 1]);
% hold on; 
% han = errorbar(dvals, dvals, dvals - d_CI2(:, 1), d_CI2(:, 2) - dvals);
% set(han, 'Color', 'b', 'LineWidth', 0.1)
% 
% xlabel('Effect size (d)');
% ylabel('Effect size (d) with CIs')
% title('Two-sample test');
% han2 = plot_horizontal_line(0); set(han2, 'Linestyle', '--');

%%



% By Tor Wager, 10/2018, alpha version. Comparisons with other methods and
% published examples should be performed, independent checks should be
% performed.

% t  = d * sqrt(N)
% we want critical value x where noncentral tcdf is 0.025 and 0.975 for 2-tailed CI
% or .05 .95 for 90% CI



if length(varargin) > 0
    dfe = varargin{1};
else
    dfe = n - 1;
end


% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------

tails = 2;
method = 'hedges';  % or 'noncentralt'
alphaval = 0.05;
n2 = NaN;           % Group 2 sample size for between-groups

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% This is a compact way to assign multiple variables. The input argument
% names and variable names must match, however:

allowable_inputs = {'tails' 'method' 'alphaval' 'n2'};

keyword_inputs = {'noncentralt'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case allowable_inputs

                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);

            case keyword_inputs
                % Skip, deal with this below

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% 2nd pass: Keyword inputs. These supersede earlier inputs
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'writetofile'
                writetofile = true;
                if isempty(movieoutfile)
                    movieoutfile = fullfile(pwd, 'rmssd_movie_file');
                    fprint('Writing file with default name:\n%s\n', movieoutfile);
                end

            case {'nomovie' 'nodisplay'}
                showmovie = false;

        end
    end
end

% Validate attributes

classes = {'numeric'};
attributes = {'scalar', 'integer', '<=',2,'>=', 1};
validateattributes(tails,classes,attributes);



% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------
switch method

    case 'hedges'
        [d_CI, ci_2sample] = effect_size_CI_hedgesolkin(d, n, n2, alphaval, tails);

        % Replace if we want the 2-sample case
        if ~isnan(n2)
            d_CI = ci_2sample; 
        end

    case 'noncentralt'

        t = d * sqrt(n);

        % Confidence interval on t
        minbnd = min(-10, 10*sign(t));
        maxbnd = max(10, 10*sign(t));

        x = [minbnd:.001:t];           % we want critical value x where noncentral tcdf is 0.025 and 0.975 for 2-tailed CI
        %p = nctcdf(t, dfe, x);  % Non-central t distribution CDF
        p = nctcdf(x, dfe, t);

        %wh = find(p < 0.975);   % 2.5% in each tail
        wh = find(p < 0.025);   % 2.5% in each tail

        t_CI = x(wh(end));    % bounds on t (noncentrality parameter)

        x = [t:.001:maxbnd];  % we want critical value x where noncentral tcdf is 0.025 and 0.975 for 2-tailed CI
        % p = nctcdf(t, dfe, x);
        p = nctcdf(x, dfe, t);

        wh = find(p < 0.975);

        t_CI(2) = x(wh(end));

        d_CI = t_CI ./ sqrt(n); % convert back to units of d


    otherwise
        error('Unknown  method')
end % case



% Get P-value and Confidence interval on P

% Convert d back to p
d2p = @(d, n, dfe, tails) tails .* (1 - tcdf(d .* sqrt(n), dfe));

P = d2p(d, n, dfe, tails);
P_CI = d2p(d_CI, n, dfe, tails);

end % main function




function [ci_1sample, ci_2sample] = effect_size_CI_hedgesolkin(d, n1, n2, alphaval, tails)
% Hedges and Olkin method
% Hedge LV, Olkin I. Statistical methods for meta-analysis. Orlando: Academic Press Inc; 2014. p. 86.
% from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5133225/
%
% Implemented here: Martin Lindquist and Tor Wager, Jan 2024
%
% One-sample example
% d = 0.5 ; n = 30;
% alphaval = 0.05;
% [ci_1sample, ci_2sample] = effect_size_CI_hedgesolkin(d, n1, n2, alphaval)

if nargin < 3  % only for stand-alone version, not needed in the subfunction
    n2 = NaN;
    alphaval = 0.05;
    tails = 2;
end

zcrit = norminv(1 - alphaval / tails);
sigma_d = sqrt(1/n1 + d.^2 / (2 * n1)); % one-sample

lb = d - zcrit * sqrt(1/n1 + d.^2 / (2 * n1));
ub = d + zcrit * sqrt(1/n1 + d.^2 / (2 * n1));

ci_1sample = [lb ub];

% Two-sample case
sigma_d = sqrt(  (n1 + n2)/(n1 * n2) + (d.^2 / (2 * (n1+n2)))  ); % two-sample

lb = d - zcrit * sigma_d;
ub = d + zcrit * sigma_d;

ci_2sample = [lb ub];

end

% NOTES from Martin - if we want to check noncentral t method further
% Computed this using two approaches:
%
% Non-centrality parameter method
%
% Lower bound:
% > pt(2.738613,29, .635) %% Non-centrality parameter obtained via trial and error to get 0.975
% [1] 0.9750502
% > .635/sqrt(30)
% [1] 0.1159346
%
% Upper bound:
% > pt(2.738613,29, 4.8) %% Non-centrality parameter obtained via trial and error to get 0.025
% [1] 0.02497157
% > 4.8/sqrt(30)
% [1] 0.8763561
%
%
% Hedges and Olkin approach:
%
% Lower bound:
% > .5-1.96*sqrt(1/30 + 0.25/60)
% [1] 0.1204476
%
% Upper bound:
% > .5+1.96*sqrt(1/30 + 0.25/60)
% [1] 0.8795524

