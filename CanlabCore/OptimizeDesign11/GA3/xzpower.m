function [varargout] = xzpower(x,varargin)
% [z80_ind,z80_grp,z80ierr,z80gerr,OUT] = xzpower(x,[c],[true],[N],[gsd],[V],[S],[trueX])
% Calculates individual and group power as a function of design matrix,
% true effect magnitudes, and contrasts
%
% Using a Monte Carlo simulation approach; xpower uses theoretical
% distribution approach
%
% tor wager, 2 / 22 / 04
%
% Defaults
% gsd = .28;          % group standard deviation, relative to true effect;
%                       default from VNL experiment (could be overestimate, bec. btwn + within
%                       not separated)
% N = 10;             % number of subjects, for group stats
% true                true effects; 1 is an increase of 1 std. dev. of scanner noise per event
%                     default is .47 SNR everywhere the first contrast has positive
%                     weights
%                     VNL estimates of signal / noise ratio for visual impulse response, 16 Hz
%                     flashing checkerboard for .5 s
%                     0.83% signal change response, 1.67% signal change error std. dev., .47 SNR 
% c = [1 1 -1 -1];    % contrast matrix; one row per contrast, one column
% per condition; default is one contrast for each effect
% 
% V                   is noise autocorrelation matrix, 
%                       symmetric and ones on the diagonal, n x n (n is number of samples)
%                       default is white noise; see also
%                       canonical_autocorrelation
%
% S                   is filtering (high/low pass) matrix; see getSmoothing
%
%
% see power_map1.m for examples and scripts
%
% outputs:
% z80_ind and z80_grp are the power estimates for individual and group
% designs.  These correspond to the z-value that is expected to be reached
% in 80% of the realizations of the design.

% Another way of creating the matrices:
%E = inline('inv(x'' * x) * x'' * Vi','x','Vi'); % E is error cov matrix
% matrix of orthogonal projection onto model space * Vi
% test: demo to verify that E is constant for constant SNR, regardless of
% sigma and t (as long as their ratio is equivalent)
% sigma = 1; t = .5; e = E(S*X,S* sigma * Vi);ee = e * e', t ./ trace(ee).^.5

% -----------------------------------------------------------------
% Set up input arguments and defaults
% -----------------------------------------------------------------

niter = 500;        % iterations of noise realizations for stats (per subject)
gsd = .28;            % group standard deviation, relative to true effect; default from VNL
N = 10;             % number of subjects, for group stats

for i = 1:length(varargin)
    switch i
        case 1
            if isempty(varargin{i})
                c = [eye(size(x,2)-1) zeros(size(x,2)-1,1)];    
                % contrast matrix; one row per contrast, one column per condition
            else
                c = varargin{i};
            end
        case 2
            if isempty(varargin{i})
                true = c(1,:) > 0;   % true effects; 1 is an increase of 1 std. dev. of scanner noise per event
                true = true .* .47;  % estimated SNR from visual experiment (VNL), n = 11
            else
                true = varargin{i};
            end
        case 3
            N = varargin{i};
        case 4
            gsd = varargin{i};
        case 5
            V = varargin{i};
        case 6
            S = varargin{i};
        case 7
            trueX = varargin{i};
    end
end

% defaults, if not entered

if ~(exist('c') == 1), c = [eye(size(x,2)-1) zeros(size(x,2)-1,1)];   ,end
if ~(exist('true') == 1),true = 0.47 .* (c(1,:) > 0);, end

if ~(exist('V') == 1), V = eye(size(x,1)); xc = [1 0 0];
    %[xc,V] = canonical_autocorrelation    
end
if ~(exist('S') == 1), S = []; ,end
    
if ~(exist('trueX') == 1) || isempty(trueX), trueX = x;, end

% defaults, if empty
if isempty(true),true = 0.47 .* (c(1,:) > 0);,end
if isempty(N),N = 10;,end
if isempty(gsd),gsd = 0.28;,end

% -----------------------------------------------------------------
% Do the power calculation - start with noncentrality parameters
% -----------------------------------------------------------------

if ~isempty(S), 
    if size(S,1) ~= size(S,2), error('S is not square.'), end
    if size(S,1) ~= size(x,1), error('S and X dims do not match.'), end
    if size(S,1) ~= size(V,1), error('S and Vi dims do not match.'), end
    
    xo = x; Vo = V;    % save original x, for true signal estimate
    x = S * x; V = S * V; 
else
    xo = x; Vo = V;
end

if size(c,2) < size(x,2), c(:,end+1) = 0;, end
if size(true,2) < size(x,2), true(:,end+1) = 0;, end

c = c';

xtxi = inv(x' * x);
sx = sqrt(diag(c' * xtxi * x' * V * V' * x * xtxi * c));

for iter = 1:niter
    
    e = mvnrnd(zeros(size(x,1),1),1); % var = 1
    %y = x * true' + V * e;  % + noisevector(size(x,1),xc,1);
    % This way, noise variance actually increases with autocorrelation, like sigma's on the diagonal
    % same as S * (yy) where yy = x * true'; but then we've filtered V, so
    % we're not adding the full noise...so we better do it the longer way
    % to be sure.
    
    % may 2006, tor replaced xo with trueX to test diff true hrfs
    y = trueX * true' + Vo * e;
    
    if ~isempty(S)
        y = S * y;  % only necessary for the "long way", otherwise we're double-filtering
    end

    [b,dev,stat] = glmfit(x(:,1:end-1),y);
    b = c' * stat.beta([2:end 1]);          % put intercept at end
    se = stat.s * sx;

    t(:,iter) = b ./ se;
    B(:,iter) = b;
    SE(:,iter) = se;
end
    
n = mean(t,2);
b = mean(B,2);
se = mean(SE,2);

% GETTING THE NONCENTRALITY PARAMETER n
% noncentrality parameter n = snr * b/se
% snr = |h| / sigma, where |h| = scaling of impulse response to neural
% event.  see Zarahn 2001, "Reference Effect for Power Analysis"
% We can simplify the computation by always assuming |h| = 1 and varying
% sigma, in the simplest case setting snr = 1.
% 
% t-values are directly related to snr and n.
% t = (|h| * n) / sigma
% in the simplest case, where snr = 1, t = n.
%
% n is contrast weight invariant.  contrasts can be normalized by sum of
% squared contrast weights.
%
% FROM n to POWER
% By the central limit theorem, n = t if snr = 1, and n is normally
% distributed with std 1.  Then we can integrate the normal distribution.
%
% but since the noncentrality param is essentially a t-value, as is the 20%
% level for the 80% power value, we have to do a small additional
% adjustment: convert t to Z scores based on the df.

%t80_ind = norminv(.2,n,1);  %"theoretical" estimate used in xpower.m

% Degrees of freedom: Worsley and Friston (1995) way (effective df), 
% df reduction due to autocorrelation.  This will be different if 
% a prewhitening strategy is used.

R = eye(size(x,1)) - x * xtxi * x';
U = R * V * V';
df = (trace(U) .^ 2) ./ (trace(U * U));

%df = size(x,1) - size(x,2);  % standard way, but not correct with
%autocorrelated errors

t80_ind = prctile(t',20)';
tmp = norminv(tcdf(t80_ind,df));
tmp = min(tmp,t80_ind);             % use min, as conversion fails for high-power
    
varargout{1} = tmp;

% GROUP POWER
% group power is calculated by considering the contribution of the
% individual se to the group se, and assuming a fixed intersubject variance (set nominally to 1).  with a
% with group rrue intersubject standard deviation (gsd) of 1, then
% gse = (gsd + E(indiv_se)) ./ sqrt(2N)
% by the central limit theorem.
% we use sqrt(2N) because std dev of n rnd variables = (s1 + s2 ... sn) / sqrt(n)
% individual se (ise) is the standard deviation of the sampling distribution
% gse is the standard deviation
% 
% take as ise the coefficient of variation (cv), sx / b, which is contrast
% scale independent.

if nargout > 1

    %gse = (gsd + (sx ./ b)) ./ sqrt(2 * N); % this is the theoretical version; we'll estimate here

    % take random set of N betas (with replacement), add noise
    % mvnrnd is used to generate noise, and it wants the variance, so we
    % use gsd^2 = sigma^2
    B = B'; clear bg, clear bs
    for iter = 1:niter*2
        
        B = getRandom(B); Btest = B(1:N,:); Btest = Btest + mvnrnd(zeros(size(Btest,1),1),gsd.^2);
        bg(iter,:) = mean(Btest);              % group effects for each contrast
        bs(iter,:) = std(Btest) ./ sqrt(N-1);  % group standard error  
        
    end
    
    tg = bg ./ bs;              % group t-value
    t80_grp = prctile(tg,20);   % group 80% power t-level
    
    ng = mean(tg);              % expected t-value (noncentrality); is this comparable to xpower???
    
    dfg = N - 1;
    tmp = norminv(tcdf(t80_grp,dfg));   % adjust to Z-score based on df
    tmp = min(tmp,t80_grp);             % use min, as conversion fails for high-power
    varargout{2} = tmp;
    
    
    if nargout > 2
        
        OUT.x = x;
        OUT.gsd = gsd;
        OUT.N = N;
        OUT.true = true;
        OUT.c = c';
        OUT.y = y;
        OUT.sx = sx;
        OUT.n = n;
        OUT.df = df;
        OUT.ng = ng;
        
        varargout{3} = OUT;
        
    end
    
end

return



        


