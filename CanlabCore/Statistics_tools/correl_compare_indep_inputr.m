function [rdiff, Z, pval, stats] = correl_compare_indep_inputr(r1, r2, N1, N2, varargin)
    % [rdiff, Z, pval, stats] = correl_compare_indep_inputr(r1, r2, N1, N2, varargin)
    %
    % Compare two Pearson's correlation values (r1 and r2) collected from independent
    % samples with sample sizes N1 and N2
    %
    % Tor Wager, March 2010
    %
    % Based on: www.stat-help.com/
    % see also: correl_compare_indep.m and correl_compare_dep.m
    %
    % Example:
    % [rdiff, Z, pval, stats] = correl_compare_indep_inputr(.33, .77, 100, 100)
    
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

