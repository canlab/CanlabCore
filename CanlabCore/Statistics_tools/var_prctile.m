function [var_prc,pci,Nneeded,pcitarget] = var_prctile(p,n_in_sample,varargin)
% :Usage:
% ::
%
%     [var_prc,pci,Nneeded,pcitarget] = var_prctile(p,n_in_sample,['nboot',nboot],['x',data])
%
%     [var_prc,pci,Nneeded,pcitarget] = var_prctile(p,50,'nboot',5000)
%
% :Inputs:
%
%   **p:**
%        is desired prctile of data (threshold)
%
%   **nboot:**
%        is number of bootstrap samples (if bootstrapping)
%
%   **n_in_sample:**
%        is number of obs. in original sample
%
%   **x:**
%        is data sample of distribution of interest (empirical pdf based on
%        this)
%
%        Note: using empirical PDF/CDF depends a great deal on choice of h (see
%        code).
%
% :Examples:
% ::
%
%    Nneeded = [];
%    for p = [.05:-.001:.001]
%        [var_prc,pci,Nneeded(end+1),pcitarget] = var_prctile(x,p);
%    end
%    figure;
%    plot([.05:-.001:.001],Nneeded)
%
% ..
%    tor wager, jan 2007
% ..

nboot = Inf;
norm_model = 1;
x = [];
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'nboot', nboot = varargin{i+1};                    % # bootstrap samples
            case {'x','data'}, x = varargin{i+1}; norm_model = 0;   % use empirical PDF

            otherwise
                error('Unknown string option');
        end
    end
end

% ------------------------------------------------------------
% * get PDF
% ------------------------------------------------------------


if norm_model
    % get normal PDF at p-th percentile
    % -------------------------------------------
    prc = norminv(p);
    pdfval = normpdf(prc);
else
    % empirical estimate of PDF from x (data)
    % -------------------------------------------
    % prc is x-score at prctile of interest
    prc = prctile(x,100*p);

    % now we need pdf at that prctile.
    % pick some unit h, and differentiate around prc

    h = max(.01,1000 ./ length(x));  % enough so we have reasonable idea
    h = min(h,p);
    h = max(h,1000 ./ length(x));
    
    pdfval = 0;

    while pdfval == 0
        pdfval = emppdf(x,prc,h);

        h = h + .01;
    end
end


% ------------------------------------------------------------
% * get variance
% ------------------------------------------------------------

%for normal: fit = ( 1./ (normpdf(norminv(p)).^2) ) .* (p *(1-p)./N);
% from Brown and Wolfe, asymptotic
%var_prc = ( 1./ (pdfval.^2) ) .* (p *(1-p)./N);

% from Martin's bootstrap book, Ch 19, Eq. 19.7, p. 275
var_prc = ( 1./ (pdfval.^2) ) .* (p *(1-p)) .* (1./nboot + 1./n_in_sample);

if nargout < 2, return, end


% ------------------------------------------------------------
% * get confidence interval
% ------------------------------------------------------------

halfci_prc = 1.96 .* sqrt(var_prc);

% lower and upper bounds on p-value derived from distribution of x
pci = get_pval_ci(x,prc,halfci_prc);

if nargout < 3, return, end

% now get nboot needed to achieve target
% ---------------------------------------

% expression for tolerance: if upper bound of p-value is within 10% of p-value
maxN = 50000;               % maximum number of iterations
Nneeded = 500;              % starting estimate for nboot needed (should be fairly low)
Nstep = 500;                % increase Nneeded in units of Nstep until satisfied
p_upper_bound = p + .1 * p; % upper bound for p desired; now 10% larger than p
pcitarget = Inf;            % target conf. interval for p-values with Nneeded boot samples

% computations we don't have to repeat.  Divide by sqrt(N) to get halfci of prc
halfci_squared_determiner = 1.96 .* sqrt( ( 1./ (pdfval.^2) ) .* p *(1-p) );

while ( max(pcitarget) > p_upper_bound ) && ( Nneeded < maxN )

    Nneeded = Nneeded + Nstep;
    halfci_prc = halfci_squared_determiner ./ sqrt(Nneeded);
    pcitarget = get_pval_ci(x,prc,halfci_prc);

end





end



function pci = get_pval_ci(x,prc,halfci_prc)
% get confidence interval for p-value, based on data x, x-axis-value prc,
% and confidence half-interval for prc

if isempty(x)
    % normal CDF
    pci = [max(0,normcdf(prc - halfci_prc)) min(1,normcdf(prc + halfci_prc))];
else
    % empirical CDF
    pci = [max(0,empcdf(x,prc - halfci_prc)) min(1,empcdf(x,prc + halfci_prc))];
end

end



function pdfval = emppdf(x,prc,h)
% empirical PDF of data numerically differentiated in a window h
pdfval = ( sum(x <= prc + h) - sum(x <= prc - h) ) ./ (length(x) * 2 * h);

end


function p = empcdf(x,prc)
% empirical CDF of data x at prctile prc

p = sum(x <= prc) ./ length(x);
end
