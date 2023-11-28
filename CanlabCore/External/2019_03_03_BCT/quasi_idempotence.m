function [XN,IOTA,EPS,U]=quasi_idempotence(X,K)
%QUASI_IDEMPOTENCE     Connection matrix quasi-idempotence
%
%   [XN,IOTA,EPS,U]=quasi_idempotence(X)
%   [XN,IOTA,EPS,U]=quasi_idempotence(X,K)
%
%   The degree of quasi-idempotence of a matrix represents how close it 
%   is to being idempotent, i.e. invariant to squaring. In turn, this 
%   reflects how closely related the edges in the graph it corresponds to 
%   are to the sums of all triangles between the corresponding nodes, 
%   spanning the entirety of the network. This probes a form of 
%   self-similarity, intended not as between nodes and modules, but 
%   between as edges and triangles (or more generally paths).
%   Networks wherein the edge strengths are weakly related to the 
%   nodal strengths have low quasi-idempotence, and networks wherein the 
%   edge strengths are significantly influenced by the strengths of the
%   corresponding nodes have high quasi-idempotence. In other words the
%   degree of quasi-idempotence may be viewed as a measure of
%   "collectivity", meaning how strongly individual nodes participate
%   in collective dynamics.
%   
%   Inputs:
%       X,
%           undirected weighted/binary connection matrix with 
%           non-negative weights,
%       K,
%           number of iterations (optional),
%               K<inf,        iterate a predetermined number of times,
%               K=inf,        iterate until convergence is attained (default).
%
%       Note: see Minati et al. (2017) for a discussion of the issues
%       associated with applying this measure to binary graphs.
%
%   Outputs:
%       XN,
%           the final matrix, from the last iteration,
%       IOTA,
%           the vector of correlation coefficients, one per iteration,
%       EPS,
%           the vector of errors, one per iteration,
%       U,
%           the number of iterations performed.
%
%       Note: see Minati et al. (2017) for a discussion of the 
%       significance of IOTA(1) and IOTA(end).
%
%   Example:
%       % Create a random undirected weighted connection matrix
%       N=100;
%       X=rand(N);
%       X=X+X';
%       [XN,IOTA,EPS,U]=quasi_idempotence(X);
%       % observe that IOTA~0
%       disp(IOTA);
%       % 
%       % make the connection matrix idempotent by squaring it
%       X=X^2;
%       [XN,IOTA,EPS,U]=quasi_idempotence(X);
%       % observe that IOTA~1
%       disp(IOTA);
%
%   Reference:
%       Minati et al. (2017) CHAOS, 043115
%
%   Ludovico Minati, University of Trento, Trento, Italy
%   and Institute of Nuclear Physics, Polish Academy of Sciences, Krakow,
%   Poland, 2017-2018
%   lminati@ieee.org
%
%   Modification history
%   2018: Original

X=cast(X,'double'); % enforce double-precision format
if any(any(isnan(X))) % bail out if any nan found
    error('found nan, giving up!');
end
if any(any(X<0)) % bail out if any negative elements found
    error('found a negative element, giving up!');
end
if ~exist('K','var')||isempty(K)
    K=inf;
end
N=length(X); % get matrix size
X(eye(N)>0)=0; % null the diagonal in the initial matrix
X=X/norm(X); % normalize to unit norm
XN=X;
mask=triu(ones(N),1)>0; % create mask for superdiagonal elements
U=0; % initialize iterations counter
IOTA=[]; % this vector will contain the correlation coefficients
EPS=inf; % and this will contain the errors
if isinf(K)
    while EPS(end)>eps('double') % iterate until error below precision
        U=U+1; % increase iteration counter
        XN_hat=XN; % save the initial matrix
        XN=XN^2; % square the matrix
        XN=XN/norm(XN); % normalize it again (for numerical reasons)
        IOTA(end+1)=corr(X(mask),XN(mask)); % calculate correlation coefficient
        EPS(end+1)=norm(XN_hat-XN); % calculate error
    end
else
    while U<K % iterate a prescribed number of times
        U=U+1; % increase iteration counter
        XN_hat=XN; % save the initial matrix
        XN=XN^2; % square the matrix
        XN=XN/norm(XN); % normalize it again (for numerical reasons)
        IOTA(end+1)=corr(X(mask),XN(mask)); % calculate correlation coefficient
        EPS(end+1)=norm(XN_hat-XN); % calculate error
    end
end
EPS(1)=[];
