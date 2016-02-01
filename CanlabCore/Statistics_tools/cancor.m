function [cc,stats] = cancor(X,w,varargin)
% :Usage:
% ::
%
%     [cc,stats] = cancor(X,w,[permutations],[MCD robust outlier removal])
%
% :Inputs:
%
%   **X:**
%        data matrix, columns are variables, rows observations
%
%   **w:**
%        first w columns are Set 1, rest are Set 2
%
%   **cc:**
%        canonical correlations
%
%   **UV:**
%        canonical variates
%
%   **UVa:**
%        caonical variates for set a
%
%   **UVb:**
%        canonical variates for set b
%
%   **ccor:**
%        correlations between canonical variates and X variables
%
%   **ab:**
%        weights of canon. variates in row vectors, e.g., U = Xa'
%
% Johnson & Wichern, Applied Multivariate Statistical Analysis, 5th ed. 
% tested on examples from the book.  cc and ab are right, not positive about ccor
%
% :Example:
%
% Compute canonical correlations between first 2 and last 2 columns
% of X, testing against permuted column data with 1000 iterations.
% ::
%
%    [cc,stats] = cancor(X,2,1000,1);
%
% ..
%    tor wager, 6 / 4 / 03
% ..

if length(varargin) > 0, perms = varargin{1};,else,perms = 0;,end

% ----------------------------------------------------------
% outliers
% ----------------------------------------------------------

if length(varargin) > 1
        [res]=fastmcd_noplot(X);
        
        % remove n most extreme outliers and recompute correlation
        wh = res.flag==0; nout = sum(res.flag==0);
        X(wh,:) = [];
        stats.nout = nout;
end

% ----------------------------------------------------------
% canonical variates
% ----------------------------------------------------------

[cc,UV,ccor,ab,UVa,UVb] = docor(X,w);

for i = 1:perms,   % permute columns and test H0: no relation
    
    xtst=X;
    for j = 1:size(X,2)
        xtst(:,j) = getRandom(X(:,j));
    end
    cctst(i,:) = docor(xtst,w);
    
end
   
stats.cc = cc;
stats.UV = UV;
stats.UVa = UVa;
stats.UVb = UVb;
stats.ab = ab;
stats.ccor = ccor;

if perms
    stats.cctst = cctst;
    stats.perms = perms;
    stats.thresh = prctile(cctst,95);
    stats.sig = cc > stats.thresh;
    stats.p = sum(cctst > repmat(cc,size(cctst,1),1)) ./ size(cctst,1);
end

return



function [cc,UV,ccor,ab,UVa,UVb] = docor(X,w)

nomore = min(w,size(X,2)-w);

c = corrcoef(X);

c11 = c(1:w,1:w);
c12 = c(1:w,w+1:end);
c21 = c(w+1:end,1:w);
c22 = c(w+1:end,w+1:end);

m1 = c11^-.5 * c12 * c22^-1 * c21 * c11^-.5;
%m2 = c22^-.5 * c21 * c11^-1 * c12 * c22^-.5;

%[e,d1]=eig(m1); % this changes the order of the eigenvectors arbitrarily, very annoying
[e,d1]=pcacov(m1);

e = e(:,1:nomore);
d1 = d1(1:nomore);
a = c11^-.5  * e;

b = c22^-1 * c21 * a;       % unscaled version of b
%tmp2 = tmp' * c22 * tmp    % = d1

%cc = sqrt(diag(d1))';        % canonical correlations
cc = sqrt(d1)';
cv=repmat(1./cc,size(b,1),1);
b =  b .* cv;

% wrong
%[f,d2]=eig(m2);
%b = c22^-.5 * f;

if nargout > 1
    
    UV = X * [a;b];

    UVa = X(:,1:length(a)) * a;
    UVb = X(:,length(a)+1:end) * b;
    
    % correlations of X with canonical
    % wrong for some reason
    V = diag(std(X)); 
    cor1 = a' * c11 * V(1:w,1:w); 
    cor2 = b' * c22 * V(w+1:end,w+1:end); 
    ccor = [cor1 cor2];
    
    for i = 1:size(UV,2)
        for j = 1:size(X,2)
            tmp = corrcoef(UV(:,i),X(:,j));
            ccor(i,j) = tmp(1,2);
        end
    end
    
    ab = [a;b]';
end

return
