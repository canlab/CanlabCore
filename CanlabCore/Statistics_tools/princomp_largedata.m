function [coeff, score, latent, tsquare] = princomp_largedata(x,econFlag)
% Principal Components Analysis
%
% This is a version created for large data sets by Matthew Davidson
%
% The default Matlab PRINCOMP is naive to large data sets. Out-of-memory
% errors are easily obtained on imaging data. The solution is to replace
% concise but inefficient calls to repmat with loops and to eliminate
% large, unused variables.
%
% Tested on random data sets and produces identical output as original
% PRINCOMP. Speed penalty is drastic for small data. 50% slower on
% 50x1000 element data set, but on a 50x10000, only 2% slower.
%
% COEFF = PRINCOMP(X) performs principal components analysis on the N-by-P
% data matrix X, and returns the principal component coefficients, also
% known as loadings.  Rows of X correspond to observations, columns to
% variables.  COEFF is a P-by-P matrix, each column containing coefficients
% for one principal component.  The columns are in order of decreasing
% component variance.
%
% PRINCOMP centers X by subtracting off column means, but does not
% rescale the columns of X.  To perform PCA with standardized variables,
% i.e., based on correlations, use PRINCOMP(ZSCORE(X)).  To perform PCA
% directly on a covariance or correlation matrix, use PCACOV.
%
% [COEFF, SCORE] = PRINCOMP(X) returns the principal component scores,
% i.e., the representation of X in the principal component space.  Rows
% of SCORE correspond to observations, columns to components.
%
% [COEFF, SCORE, LATENT] = PRINCOMP(X) returns the principal component
% variances, i.e., the eigenvalues of the covariance matrix of X, in
% LATENT.
%
% [COEFF, SCORE, LATENT, TSQUARED] = PRINCOMP(X) returns Hotelling's
% T-squared statistic for each observation in X.
%
% When N <= P, SCORE(:,N:P) and LATENT(N:P) are necessarily zero, and the
% columns of COEFF(:,N:P) define directions that are orthogonal to X.
%
% [...] = PRINCOMP(X,'econ') returns only the elements of LATENT that are
% not necessarily zero, i.e., when N <= P, only the first N-1, and the
% corresponding columns of COEFF and SCORE.  This can be significantly
% faster when P >> N.
%
% :See Also: BARTTEST, BIPLOT, CANONCORR, FACTORAN, PCACOV, PCARES, ROTATEFACTORS.
%
% :References:
%   1. Jackson, J.E., A User's Guide to Principal Components,
%      Wiley, 1988.
%   2. Jolliffe, I.T. Principal Component Analysis, 2nd ed.,
%      Springer, 2002.
%   3. Krzanowski, W.J., Principles of Multivariate Analysis,
%      Oxford University Press, 1988.
%   4. Seber, G.A.F., Multivariate Observations, Wiley, 1984.
%
% ..
%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 2.9.2.9 $  $Date: 2006/10/02 16:35:01 $
% ..

% ..
%    When X has more variables than observations, the default behavior is to
%    return all the pc's, even those that have zero variance.  When econFlag
%    is 'econ', those will not be returned.
% ..

if nargin < 2, econFlag = 0; end

[n,p] = size(x);
if isempty(x)
    pOrZero = ~isequal(econFlag, 'econ') * p;
    coeff = zeros(p,pOrZero); coeff(1:p+1:end) = 1;
    score = zeros(n,pOrZero);
    latent = zeros(pOrZero,1);
    tsquare = zeros(n,1);
    return
end

% Center X by subtracting off column means
%x0 = x - repmat(mean(x,1),n,1);
for i=1:size(x, 2)
    x(:,i) = x(:,i) - mean(x(:,i));
end
r = min(n-1,p); % max possible rank of X0

% The principal component coefficients are the eigenvectors of
% S = X0'*X0./(n-1), but computed using SVD.
[U,sigma,coeff] = svd(x,econFlag); % put in 1/sqrt(n-1) later

if nargout < 2
    % When econFlag is 'econ', only (n-1) components should be returned.
    % See comment below.
    if (n <= p) && isequal(econFlag, 'econ')
        coeff(:,n) = [];
    end

else
    % Project X0 onto the principal component axes to get the scores.
    if n == 1 % sigma might have only 1 row
        sigma = sigma(1);
    else
        sigma = diag(sigma);
    end
    score = U .* repmat(sigma',n,1); % == x*coeff
    sigma = sigma ./ sqrt(n-1);

    % When X has at least as many variables as observations, eigenvalues
    % n:p of S are exactly zero.
    if n <= p
        % When econFlag is 'econ', nothing corresponding to the zero
        % eigenvalues should be returned.  svd(,'econ') won't have
        % returned anything corresponding to components (n+1):p, so we
        % just have to cut off the n-th component.
        if isequal(econFlag, 'econ')
            sigma(n,:) = []; % make sure this shrinks as a column
            coeff(:,n) = [];
            score(:,n) = [];

        % Otherwise, set those eigenvalues and the corresponding scores to
        % exactly zero.  svd(,0) won't have returned columns of U
        % corresponding to components (n+1):p, need to fill those out.
        else
            sigma(n:p,1) = 0; % make sure this extends as a column
            score(:,n:p) = 0;
        end
    end

    % The variances of the pc's are the eigenvalues of S = X0'*X0./(n-1).
    latent = sigma.^2;

    % Hotelling's T-squared statistic is the sum of squares of the
    % standardized scores, i.e., Mahalanobis distances.  When X appears to
    % have column rank < r, ignore components that are orthogonal to the
    % data.
    if nargout == 4
        if n > 1
            q = sum(sigma > max(n,p).*eps(sigma(1)));
            if q < r
                warning('stats:princomp:colRankDefX', ...
                        ['Columns of X are linearly dependent to within machine precision.\n' ...
                         'Using only the first %d components to compute TSQUARED.'],q);
            end
        else
            q = 0;
        end
        tsquare = (n-1) .* sum(U(:,1:q).^2,2); % == sum((score*diag(1./sigma)).^2,2)
    end
end
