function [rci,sig,Z,p,rcrit] = r2z(r,n,varargin)
% Fisher's r to Z transformation for providing CIs for a correlation
%
% :Usage:
% ::
%
%     [rci,sig,Z,p,rcrit] = r2z(r,n,[alph])
%
% :Inputs:
%
%   **n:**
%        n, # of observations going into correlation
%        (count each row/subject once)
%
%        df = n - 3
%
%   **alph:**
%        two-tailed p-value cutoff, 
%        default is p < .05 
%
% :Outputs:
%
%   **rci:**
%        confidence interval in correlation values
%
%   **sig:**
%        significant at alpha value?
%
%   **Z:**
%        z-scores of correlations
%
%   **p:**
%        p-values
%
% can take a vector of r values
%
% :Examples:
% ::
%
%    [rci,sig,z] = r2z(.1:.05:.9,5,.05); figure('Color','w');hold on; plot(.1:.05:.9,rci(:,1),'g','LineWidth',2)
%    [rci,sig,z] = r2z(.1:.05:.9,10,.05); hold on; plot(.1:.05:.9,rci(:,1),'r','LineWidth',2)
%    [rci,sig,z] = r2z(.1:.05:.9,20,.05); hold on; plot(.1:.05:.9,rci(:,1),'b','LineWidth',2)
%    [rci,sig,z] = r2z(.1:.05:.9,40,.05); hold on; plot(.1:.05:.9,rci(:,1),'m','LineWidth',2)
%    [rci,sig,z] = r2z(.1:.05:.9,80,.05); hold on; plot(.1:.05:.9,rci(:,1),'k','LineWidth',2)
%    set(gca,'FontSize',18)
%    legend({'n = 5' 'n = 10' 'n = 20' 'n = 40' 'n = 80'})
%    title('.05 Confidence interval lower bound on Pearson''s r')
%    xlabel('Correlation (r)')
%    ylabel('CI Lower Bound (r)')
%    c=.2:.01:.5;,[rci,sig,z,p,rcrit]=r2z(c,39,.05);[c' sig p]
%
%    ind=1;for i=1:10:5000,[rci,s,z,p,rc(ind)]=r2z(.5,39,.05/i);,ind=ind+1;,end
%    figure;plot(1:10:5000,rc);title('Critical r with Bonf correction'),xlabel('Comparisons')
%
% ..
%    tor wager
% ..

alph = .05;
if length(varargin) > 0, alph = varargin{1};, end

if size(r,2) > size(r,1), r = r';, end

z = .5 * log( (1+r) ./ (1-r) );

% critical Z at alpha (two-tailed)
zc = norminv(1-alph/2);

% standard deviation of z (= standard error)
sd = (1 ./ (n - 3)) .^.5;

% p-values
Z = z ./ sd; p  = (1 - normcdf(abs(Z))) .* 2;

% critical r
rcrit = (exp(2.*zc.*sd) - 1) ./ (exp(2.*zc.*sd) + 1);

% lower and upper bounds = confindence interval
zci(:,1) = z - (sd*zc);
zci(:,2) = z + (sd*zc);

% transform CIs back to r

rci = (exp(2.*zci) - 1) ./ (exp(2.*zci) + 1);

% cap at 1 - shoudn't need this
%rci(:,2) = min([rci(:,2) ones(size(rci,1),1)]')';

% significant?

sig = rci(:,1) > 0 | rci(:,2) < 0;

return
