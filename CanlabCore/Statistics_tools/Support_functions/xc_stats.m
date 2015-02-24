function [mxc,mlat,OUT] = xc_stats(xc,xl,varargin)
% [mxc,mlat,OUT] = xc_stats(xc,xl,varargin)
% calculate t-statistic values for a k x k x n 3-D matrix of 
% correlation (xc) and latency (xl) values across n subjects.
%
% tor wager


names = [];
if length(varargin) > 0, names = varargin{1};,end

 if isempty(xl), mlat = [];, end

% ------------------------------------------------------------
% stats on correlation values
% ------------------------------------------------------------

warning off, clear stelatency, clear tlat, clear stecorr, clear tcorr
% standard errors and t-values, element by element

% only do this on the upper triangle, to be faster

for j = 1:size(xc,2)-1
    for k = j + 1 : size(xc,2)
        
        if ~isempty(xl)
            tmp = squeeze(xl(j,k,:));
            z = tmp;     
            [stelatency(j,k),tlat(j,k)] = ste(z);   % t-test on latencies; maybe should be poisson dist?
        end
        
        tmp = squeeze(xc(j,k,:));
        z = .5 * log( (1+tmp) ./ (1-tmp) );     % Fisher's r-to-z transform
        [stecorr(j,k),tcorr(j,k)] = ste(z);     % correl in rand fx analysis across ss
    end
end
warning on


% makes symmetric matrices from upper tri

if ~isempty(xl), stelatency(end+1,:) = 0; ,end
stecorr(end+1,:) = 0;
tcorr(end+1,:) = 0; 
if ~isempty(xl), 
    tlat(end+1,:) = 0; 
    stelatency = stelatency + stelatency';
end

stecorr = stecorr + stecorr' + Inf * eye(size(stecorr,1));
tcorr = tcorr + tcorr';     % t-scores for significance/reliability of cross-correlations across subjects

if ~isempty(xl)
    tlat = tlat + tlat';        % same for time lags across subjects
end

% group means

if ~isempty(xl)
    mlat = mean(xl,3);          % mean of latency data (cross-lag, group average)
end
mxc = mean(xc,3);           % mean cross-correlations among k regions, group average


% ------------------------------------------------------------
% print and save output
% ------------------------------------------------------------


OUT.desc1 = 'xl is latency mtx for each subject, xc is individual correlation mtx';
if ~isempty(xl), OUT.xl = xl;, end
OUT.xc = xc;
OUT.desc2 = 'mlat and mxc: matrices of means across subjects for latency and correl'
OUT.desc3 = 'tlat/tcorr are t-values across z-transformed correlations';
OUT.mxc = mxc;

if ~isempty(xl)
    OUT.mlat = mlat;
    OUT.tlat = tlat;
    OUT.tcorr = tcorr;
end

% Output
if ~isempty(xl)
    fprintf(1,'Latency\n');
    str = correlation_to_text(mlat,Inf,names);
end

fprintf(1,'correlation T-values, p<.05 corrected, assuming normal distribution of latency diffs\n');
% Alpha correction - bonferroni.
numobs=(size(mxc,1)*(size(mxc,1)-1))/2;
corrp=1-(0.05/ (2 * (   numobs   )));       % 2-tailed corr p
crit_t = tinv_t(corrp,size(xc,1)-3);          % critical t-value for corrected significance
crit_tu = tinv_t(1-(.05/2),size(xc,1)-3);         % critical t-value for uncorrected significance

if ~isempty(xl)
    str = correlation_to_text(tlat,crit_t,names);
end

OUT.crit_t = crit_t;
OUT.crit_tu = crit_tu;
OUT.crit_tdesc = 'crit_t crit t value corrected, _tu uncorrected';
OUT.corrp = corrp;
OUT.corrpdesc = 'critical corrected p-value (bonferroni)';


fprintf(1,'Average correlation across subjects\n');
[str, sigu] = correlation_to_text(mxc,Inf,names);

fprintf(1,'Corr T-values, p<.05 uncorrected, 2-tailed, assuming independent samples over time\n');
[str, sigu] = correlation_to_text(tcorr,crit_tu,names);
OUT.sigmat_uncorrected = sigu;

fprintf(1,'Corr T-values, p<.05 corrected, 2-tailed, assuming independent samples over time\n');
[str, sigmat] = correlation_to_text(tcorr,crit_t,names);

OUT.sigmat = sigmat;




return