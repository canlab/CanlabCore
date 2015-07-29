
function a = bagging(alg,hyper) 

% ============================================================================
% BAGGING bagging object
% ============================================================================
% A=bagging(A,H) returns a bagging object initialized with algorithm A and
%                hyperparameters H.
%
% Bagging trains BAGS classifiers on resamplings with replacement of the data
% of size M.  The final classifier is the averaged classifier.
% For classifiers that output a real-valued output in pattern recognition this
% gives you the choice of performing the average before or after taking the
% sign.
%
% You can also use the bagging object to average the results of already
% trained classifiers, simply define BAGGING(A) where A is a GROUP of
% classifiers 
%
% Hyperparameters (with defaults)
%   child=svm            -- the algorithm to bag
%   bags=10              -- the number of times to bag
%   m=500                -- the number of points to sample for each bag
%                           (can also be a fraction of the training data)
%
% Example:
%    a1=svm; a1.C=10;
%    a=bagging(a1); a.bags=10; a.m=20;
%    [r a]=train(a,toy2d);
%    loss(test(a,toy2d))
%
% ============================================================================
% Reference : Bagging Predictors
% Author    : Leo Breiman
% Link      : http://citeseer.lcs.mit.edu/breiman96bagging.html
% ============================================================================

a.bags=10; a.m=500; a.child=svm;

p=algorithm('bagging');
a=class(a,'bagging',p);

if nargin>0
	a.child=alg;
end

if strcmp(a.child.algorithm.name,'group')
	%% already defined machines, can choose bags automatically
	a.bags=length(a.child.child);
end

if nargin==2
	eval_hyper;
end

