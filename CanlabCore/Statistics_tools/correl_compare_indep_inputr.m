function [rdiff, Z, pval, stats] = correl_compare_indep_inputr(r1, r2, N1, N2, varargin)
% Compare two Pearson's correlation values collected from independent samples
% :Based on: www.stat-help.com/
% :Usage:
% ::
%     [rdiff, Z, pval, stats] = correl_compare_indep_inputr(r1, r2, N1, N2, varargin)
%
% :Inputs:
%
%   **r1:**
%        Correlation coefficient for group 1
%
%   **r2:**
%        Correlation coefficient for group 2
%
%   **N1:**
%        Sample size for group 1  
%
%   **N2:**
%        Sample size for group 2 
%
% :Outputs:
%
%   **rdiff:**
%        Difference in correlation coefficient matrices (r1 - r2)
%
%   **Z:**
%        Z score for difference in correlation coefficients.
%
%   **pval:**
%        p value for two-tailed z-test
%
%   **stats:**
%        Structure of variables including r1, r2, N1, N2, rdiff, zdiff, Z,
%        pval, myalpha, sig
%
% :Examples:
% ::
%
%    %generate two random input matrices (20 subjects x 5 variables)
%    r1 = .77;
%    r2 = .33;
%    N1 = 100;
%    N2 = 100;
%
%   [rdiff, Z, pval, stats] = correl_compare_indep_inputr(r1, r2, N1, N2);
%
%
% :See also:
%   - correl_compare_indep_inputr, correl_compare_permute*
%
%
% ..
%    Tor Wager, March 2010
% ..

    if nargin == 0
        disp('Using sample values for r1 r2 N1 N2')
        N1 = 100;
        N2 = 100;
        r1 = .33;
        r2 = .77;
    end

    rdiff = r1 - r2;

    sediff = sqrt(1/(N1-3) + 1/(N2-3));  % Standard error of difference between two standard normals, with df = N - 3 for correlation coeff
    zdiff = fisherz(r1) - fisherz(r2);  % difference in Fisher's z values

    Z = zdiff ./ sediff;   % Z-score (inferential stat)

    pval = 2 * (1 - normcdf(abs(Z))); % p-value

    myalpha = .05;  % ROI
    sig = pval < myalpha;

    stats = struct('r1', r1, 'r2', r2, 'N1', N1, 'N2', N2, 'diff_r_values', rdiff, ...
        'fishers_z_diff', zdiff, 'Z_stat', Z, 'p', pval, 'alpha', myalpha, 'sig', sig);

end

%out = correl_compare_dep(y1,y2, 'alpha',myalpha)

