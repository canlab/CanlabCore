function Y = pdist(X,s,varargin)
%PDIST Pairwise distance between observations.
%   Y = PDIST(X) returns a vector Y containing the Euclidean distances
%   between each pair of observations in the M-by-N data matrix X.  Rows of
%   X correspond to observations, columns correspond to variables.  Y is an
%   (M*(M-1)/2)-by-1 vector, corresponding to the M*(M-1)/2 pairs of
%   observations in X.
%
%   Y = PDIST(X, DISTANCE) computes Y using DISTANCE.  Choices are:
%
%       'euclidean'   - Euclidean distance
%       'seuclidean'  - Standardized Euclidean distance, each coordinate
%                       in the sum of squares is inverse weighted by the
%                       sample variance of that coordinate
%       'cityblock'   - City Block distance
%       'mahalanobis' - Mahalanobis distance
%       'minkowski'   - Minkowski distance with exponent 2
%       'cosine'      - One minus the cosine of the included angle
%                       between observations (treated as vectors)
%       'correlation' - One minus the sample correlation between
%                       observatons (treated as sequences of values).
%       'hamming'     - Hamming distance, percentage of coordinates
%                       that differ
%       'jaccard'     - One minus the Jaccard coefficient, the
%                       percentage of nonzero coordinates that differ
%       function      - A distance function specified using @, for
%                       example @DISTFUN
%
%   A distance function must be of the form
%
%         function D = DISTFUN(XI, XJ, P1, P2, ...),
%
%   taking as arguments two L-by-N matrices XI and XJ each of which
%   contains rows of X, plus zero or more additional problem-dependent
%   arguments P1, P2, ..., and returning an L-by-1 vector of distances D,
%   whose Kth element is the distance between the observations XI(K,:)
%   and XJ(K,:).
%
%   Y = PDIST(X, DISTFUN, P1, P2, ...) passes the arguments P1, P2, ...
%   directly to the function DISTFUN.
%
%   Y = PDIST(X, 'minkowski', P) computes Minkowski distance using the
%   positive scalar exponent P.
%
%   The output Y is arranged in the order of ((1,2),(1,3),..., (1,M),
%   (2,3),...(2,M),.....(M-1,M)), i.e. the upper right triangle of the full
%   M-by-M distance matrix.  To get the distance between the Ith and Jth
%   observations (I < J), either use the formula Y((I-1)*(M-I/2)+J-I), or
%   use the helper function Z = SQUAREFORM(Y), which returns an M-by-M
%   square symmetric matrix, with the (I,J) entry equal to distance between
%   observation I and observation J.
%
%   Example:
%
%      X = randn(100, 5);                 % some random points
%      Y = pdist(X, 'euclidean');         % unweighted distance
%      Wgts = [.1 .3 .3 .2 .1];           % coordinate weights
%      Ywgt = pdist(X, @weucldist, Wgts); % weighted distance
%
%      function d = weucldist(XI, XJ, W) % weighted euclidean distance
%      d = sqrt((XI-XJ).^2 * W');
%
%   See also SQUAREFORM, LINKAGE, SILHOUETTE.

%   An example of distance for data with missing elements:
%
%      X = randn(100, 5);     % some random points
%      X(unidrnd(prod(size(X)),1,20)) = NaN; % scatter in some NaNs
%      D = pdist(X, @naneucdist);
%
%      function d = naneucdist(XI, XJ) % euclidean distance, ignoring NaNs
%      sqdx = (XI-XJ).^2;
%      nk = sum(~isnan(sqdx),2); nk(nk == 0) = NaN;
%      d = sqrt(nansum(sqdx')' .* size(XI,2) ./ nk); %correct for missing coords

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 1.15 $

if nargin < 2
    s = 'euc';
    distfun = @distcalc;
    distargs = {s};
else
    if ischar(s)
        methods = strvcat('euclidean','seuclidean','cityblock','mahalanobis','minkowski','cosine','correlation','hamming','jaccard');
        i = strmatch(lower(s), methods);
        if length(i) > 1
            error(sprintf('Ambiguous ''DISTANCE'' argument:  %s.', s));
        elseif isempty(i)
            % error(sprintf('Unknown ''DISTANCE'' argument:  %s.', s));
            distfun = str2func(s);
            distargs = varargin;
            s = 'usr';
        else
            s = lower(methods(i,1:3));
            distfun = @distcalc;
            distargs = {s};
        end
    elseif isa(s, 'function_handle') |  isa(s, 'inline')
        distfun = s;
        distargs = varargin;
        s = 'usr';
    else
        error('The ''DISTANCE'' argument must be a string or a function.');
    end
end

[m, n] = size(X);
if any(imag(X(:))) & ~isequal(s,'usr')
   error('PDIST does not accept complex inputs.');
end

switch s
case 'seu' % Standardized Euclidean weights by coordinate variance
   distargs{end+1} = 1 ./ var(X)';
case 'mah' % Mahalanobis
   distargs{end+1} = inv(cov(X));
case 'min' % Minkowski distance needs a third argument
   if nargin < 3  % use default value for exponent
      distargs{end+1} = 2;
   elseif varargin{1} > 0
      distargs{end+1} = varargin{1}; % get exponent from input args
   else
      error('The exponent for the Minkowski metric must be positive.');
   end
case 'cos' % Cosine
   Xnorm = sqrt(sum(X.^2, 2));
   if min(Xnorm) <= eps * max(Xnorm)
       error(['Some points have small relative magnitudes, making them ', ...
              'effectively zero.\nEither remove those points, or choose a ', ...
              'distance other than cosine.'], []);
   end
   X = X ./ Xnorm(:,ones(1,n));
case 'cor' % Correlation
   X = X - repmat(mean(X,2),1,n);
   Xnorm = sqrt(sum(X.^2, 2));
   if min(Xnorm) <= eps * max(Xnorm)
       error(['Some points have small relative standard deviations, making ', ...
              'them effectively constant.\nEither remove those points, or ', ...
              'choose a distance other than correlation.'], []);
   end
   X = X ./ Xnorm(:,ones(1,n));
end

if m < 2
   % Degenerate case, just return an empty of the proper size.
   Y = zeros(1,0);
   return;
end

% Create (I,J) defining all pairs of points
p = (m-1):-1:2;
I = zeros(m*(m-1)/2,1);
I(cumsum([1 p])) = 1;
I = cumsum(I);
J = ones(m*(m-1)/2,1);
J(cumsum(p)+1) = 2-p;
J(1)=2;
J = cumsum(J);

% For large matrices, process blocks of rows as a group
n = length(I);
ncols = size(X,2);
blocksize = 1e4;                     % # of doubles to process as a group
M = max(1,ceil(blocksize/ncols));    % # of rows to process as a group
nrem = rem(n,M);
if nrem==0, nrem = min(M,n); end

Y = zeros(1,n);
ii = 1:nrem;
try
    Y(ii) = feval(distfun,X(I(ii),:),X(J(ii),:),distargs{:})';
catch
    if isa(distfun, 'inline')
        error(['The inline distance function generated the following ', ...
               'error:\n%s'], lasterr);
    elseif strfind(lasterr, ...
                   sprintf('Undefined function ''%s''', func2str(distfun)))
        error('The distance function ''%s'' was not found.', func2str(distfun));
    else
        error(['The distance function ''%s'' generated the following ', ...
               'error:\n%s'], func2str(distfun),lasterr);
    end
end;
for j=nrem+1:M:n
    ii = j:j+M-1;
    try
        Y(ii) = feval(distfun,X(I(ii),:),X(J(ii),:),distargs{:})';
    catch
        if isa(distfun, 'inline')
            error(['The inline distance function generated the following ', ...
                    'error:\n%s'], lasterr);
        else
            error(['The distance function ''%s'' generated the following', ...
                    'error:\n%s'], func2str(distfun),lasterr);
        end;
    end;
end

% ----------------------------------------------
function d = distcalc(XI,XJ,s,arg)
%DISTCALC Perform distance calculation for PDIST.
switch s
case 'euc',   d = sqrt(sum((XI-XJ).^2,2));            % Euclidean
case 'seu',   d = sqrt(((XI-XJ).^2) * arg);           % Standardized Euclidean
case 'cit',   d = sum(abs((XI-XJ)),2);                % City Block
case 'mah',   Y = XI - XJ;
              d = sqrt(sum((Y*arg).*Y,2));            % Mahalanobis
case 'min',   d = sum(abs((XI-XJ)).^arg,2).^(1/arg);  % Minkowski
case 'cos',   d = 1 - sum(XI.*XJ,2);                  % Cosine
case 'cor',   d = 1 - sum(XI.*XJ,2);                  % Correlation
case 'ham',   d = sum(XI ~= XJ,2) / size(XI,2);       % Hamming
case 'jac',   nz = XI ~= 0 | XJ ~= 0;
              ne = XI ~= XJ;
              d = sum(ne&nz,2) ./ sum(nz,2);          % Jaccard
end