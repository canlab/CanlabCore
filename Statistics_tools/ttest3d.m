function [mxc,t,sig,out] = ttest3d(xc)
    % [mean,t,sig,out] = ttest3d(xc)
    %
    % calculate t-statistic values for a k x k x n 3-D matrix of
    % values across n subjects.
    % stats on each element, across 3rd dim
    %
    % tor wager
    %
    % to print output, see also:
    % correlation_to_text(mxc,sig);

    % ------------------------------------------------------------
    % stats on each element, across 3rd dim
    % ------------------------------------------------------------

    warning off, clear stelatency, clear tlat, clear stecorr, clear tcorr
    % standard errors and t-values, element by element

    for j = 1 : size(xc,2)
        for k = 1 : size(xc,2)

            tmp = squeeze(xc(j,k,:));
            %z = .5 * log( (1+tmp) ./ (1-tmp) );     % Fisher's r-to-z
            %transform (old, for
            %correlation inputs)
            [stecorr(j,k), t(j,k), nobs(j, k), p(j, k), mxc(j, k)] = ste(tmp);     % correl in rand fx analysis across ss
        end
    end
    warning on

    % done above
    %mxc = nanmean(xc,3);           % mean cross-correlations among k regions, group average

    out.mean = mxc;
    out.ste = stecorr;
    out.t = t;

    out.nobs = nobs;
    
    p(p == 0) = 1000*eps; % fix for very very sig effects
    
    out.pvals = p;

    % Alpha correction - bonferroni.
    numtests = (size(mxc,1)*(size(mxc,1)-1))/2;       % number of obs in upper triangle
    corrp = 1 - (0.05 / (2 * (   numtests   )));         % 2-tailed corr p
    
    %crit_t = tinv_t(corrp,size(xc,3)-1);          % critical t-value for corrected significance
    %crit_tu = tinv_t(1-(.05/2),size(xc,3)-1);     % critical t-value for uncorrected significance

    %out.numtests = numtests;
    out.bonf_corr_pthresh = corrp;
    %out.crit_t = crit_t;
    %out.crit_tu = crit_tu;
    %out.crit_t_descrip = 'Critical t-values for Bonferoni correction on number of upper triangular elements';
    %out.crit_tu_descrip = 'Uncorrected critical t-value, p < .05, 2-tailed';
    out.numcomps = numtests;

    sig = p < corrp; %abs(t) > crit_t;
    sigu = p < .05; %abs(t) > crit_tu;

    sig = sig .* sign(mxc);
    sigu = sigu .* sign(mxc);

    sig(isnan(sig)) = 0;
    sigu(isnan(sigu)) = 0;

    out.sig = sig;
    out.sigu = sigu;

    N = size(xc, 3);
    %out.pvals = 1 - tcdf(abs(t), N - 1);
    %out.pvals_descrip = 'one-tailed p-values';

    [out.fdrthr, out.fdrsig] = fdr_correct_pvals(out.pvals, t);

end


function [pthr,sig] = fdr_correct_pvals(p,r)

    psq = p; psq(find(eye(size(p,1)))) = 0;
    psq = squareform(psq);
    pthr = FDR(p,.05);
    if isempty(pthr), pthr = 0; end

    sig = sign(r) .* (p < pthr);

    sig(isnan(sig)) = 0;
end

