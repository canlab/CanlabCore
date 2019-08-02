function out = ttest3d(xc)
% calculate within-replicate (e.g., subject) correlations and t-statistic values for 
%
% :Usage:
% ::
%
%     OUT = ttest3d(xc)
%
%
% :Inputs:
%
%   **xc:**
%        a k x k x n 3-D matrix of values across n subjects.
%
% :Outputs:
%
%   **OUT:**
%      A structure containing the following fields:
% 
%     'r:                  Correlation or mean correlation, or other association metric      '
%     'ste:                Std. error of each element of r                                   '
%     't:                  t-value across replicates (e.g., subjects) for each element of r  '
%     'nobs:               number of replicates (e.g., subjects) for each element of r       '
%     'p:                  P-value across replicates (e.g., subjects) for each element of r  '
%     'p_thr:              P-value threshold, uncorrected'
%     'sigu:               Logical matrix of significant associations, uncorrected p < 0.05  '
%     'bonf_corr_pthresh:  P-value threshold for P < 0.05 Bonferroni-corrected               '
%     'sig:                Logical matrix of significant associations, Bonferroni-corrected  '
%     'fdrthr:             P-value threshold for q < 0.05 FDR-corrected                      '
%     'fdrsig:             Logical matrix of significant associations, q < 0.05 FDR-corrected'
%     'numcomps:           Number of unique pairwise comparisons (e.g., in lower triangle)   '
%
% :Examples:
% ::
%
% Below is a complete example generating simulated data and plotting output:
% -------------------------------------------------------------------------
%
%  % simulate data matrix for 10 subjects x 20 observations each, 4 variables
% n = 10;
% covmtx = eye(4); covmtx(1, 2) = .7; covmtx(3, 4) = .7; covmtx(1, 3) = .5; covmtx(1, 4) = .5;
% covmtx = (covmtx + covmtx') ./ 2;
% x = mvnrnd([0 0 0 0], covmtx, 20*n);
% 
%  % Create subject index vector indicating subject (replicate) ID:
% s = []; for i = 1:n, s = [s; i * ones(20, 1)]; end
% 
% % % Convert to cell to get correlations within each subject:
% x_cell = canlab_mat2cell(x, s);
% 
%  % Calculate within-person correlations:
% r_cell = cellfun(@(x) corr(x, 'Type', 'Spearman'), x_cell, 'UniformOutput', false);
% r3d = cat(3, r_cell{:});
% 
%  % Do the 3-d t-test to get stats:
% OUT = ttest3d(r3d);
%
%  % Print the correlation matrix:
% correlation_to_text(OUT.r, OUT.sig);
%
%  % Plot the correlation matrix:
% plot_correlation_matrix(OUT);
%

% Programmers' notes:
% ..
%    tor wager. updated july 2019 (documentation and output format)
% ..

warning off, clear stelatency, clear tlat, clear stecorr, clear tcorr
% standard errors and t-values, element by element

p_thr = 0.05;

descrip =     {'r:                  Correlation or mean correlation, or other association metric'};
descrip(2) =  {'ste:                Std. error of each element of r'};
descrip(3) =  {'t:                  t-value across replicates (e.g., subjects) for each element of r'};
descrip(4) =  {'nobs:               number of replicates (e.g., subjects) for each element of r'};
descrip(5) =  {'p:                  P-value across replicates (e.g., subjects) for each element of r'};
descrip(6) =  {'p_thr:              P-value threshold, uncorrected'};
descrip(7) =  {'sigu:               Logical matrix of significant associations, uncorrected p < 0.05'};
descrip(8) =  {'bonf_corr_pthresh:  P-value threshold for P < 0.05 Bonferroni-corrected'};
descrip(9) =  {'sig:                Logical matrix of significant associations, Bonferroni-corrected'};
descrip(10) = {'fdrthr:             P-value threshold for q < 0.05 FDR-corrected'};
descrip(11) = {'fdrsig:             Logical matrix of significant associations, q < 0.05 FDR-corrected'};
descrip(12) = {'numcomps:           Number of unique pairwise comparisons (e.g., in lower triangle)'};
descrip = char(descrip{:});

% ------------------------------------------------------------
% stats on each element, across 3rd dim
% ------------------------------------------------------------


for j = 1 : size(xc, 1)
    for k = 1 : size(xc, 2)
        
        tmp = squeeze(xc(j, k, :));
        %z = .5 * log( (1+tmp) ./ (1-tmp) );     % Fisher's r-to-z
        %transform (old, for
        %correlation inputs)
        [stecorr(j, k), t(j, k), nobs(j, k), p(j, k), mxc(j, k)] = ste(tmp);     % correl in rand fx analysis across ss
    end
end
warning on

% done above
%mxc = nanmean(xc,3);           % mean cross-correlations among k regions, group average

out.descrip = descrip;

out.r = mxc;
out.ste = stecorr;
out.t = t;

out.nobs = nobs;

p(p == 0) = 100*eps;            % fix for very very sig effects

out.p = p;

% Alpha correction - bonferroni.
numtests = (size(mxc, 1)*(size(mxc, 2) - 1)) / 2;       % number of obs in upper triangle
corrp = 0.05 / (2 * (   numtests   ));                  % 2-tailed corr p

%crit_t = tinv_t(corrp,size(xc,3)-1);          % critical t-value for corrected significance
%crit_tu = tinv_t(1-(.05/2),size(xc,3)-1);     % critical t-value for uncorrected significance

%out.numtests = numtests;
out.bonf_corr_pthresh = corrp;
%out.crit_t = crit_t;
%out.crit_tu = crit_tu;
%out.crit_t_descrip = 'Critical t-values for Bonferoni correction on number of upper triangular elements';
%out.crit_tu_descrip = 'Uncorrected critical t-value, p < .05, 2-tailed';
out.numcomps = numtests;

sig = p < corrp;  % abs(t) > crit_t;
sigu = p < p_thr; % abs(t) > crit_tu;

sig = sig .* sign(mxc);
sigu = sigu .* sign(mxc);

sig(isnan(sig)) = 0;
sigu(isnan(sigu)) = 0;

out.sig = sig;
out.sigu = sigu;

N = size(xc, 3);
%out.p = 1 - tcdf(abs(t), N - 1);
%out.p_descrip = 'one-tailed p-values';

[out.fdrthr, out.fdrsig] = fdr_correct_pvals(out.p, t);

end


function [pthr,sig] = fdr_correct_pvals(p, t)

issquarematrix = false;

[m, n] = size(p);
if m == n, issquarematrix = true; end

psq = p;

if issquarematrix
    
    %     psq(find(eye(size(p, 1)))) = 0;
    %     psq = squareform(psq);
    
    trilp = tril(p, -1);
    wh = logical(tril(ones(size(p)), -1)); % Select off-diagonal values (lower triangle)
    psq = double(trilp(wh));             % Vectorize and enforce double
    
end

psq(psq < 10*eps) = 10*eps;        % Avoid exactly zero values

pthr = FDR(psq, .05);
if isempty(pthr), pthr = -Inf; end

sig = sign(t) .* (p <= pthr);
% use equal to because in some cases FDR will return all very low P-values,
% and all are exactly at threshold.

sig(isnan(sig)) = 0;

end

