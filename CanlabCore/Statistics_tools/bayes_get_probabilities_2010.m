function [priors_Pt, pa1_given_t, pa0_given_t, varargout] = bayes_get_probabilities_2010(Y, Xi, k)
% :Usage:
% ::
%
%     [priors_Pt, pa1_given_t, pa0_given_t, pt_given_act1, pt_given_act0, pa1_given_not_t] = bayes_get_probabilities_2010(Y, Xi, k)
%
% :Inputs:
%
%   **k:**
%        is regularization param, biases P(activity | task class) towards 0.5
%
%   **kY:**
%        is data matrix, obs x features, 1/0 (active/not)
%
%   **kXi:**
%        is task indicator matrix
%
% This is a sub-function of classify_naive_bayes.m
% For complete help, see classify_naive_bayes.m
%
% ..
%    tor wager, oct 07
% ..

[nobs, nfeatures] = size(Y);
nclasses = size(Xi, 2);

% priors_Pt and evidence: Marginal probability estimates
% -----------------------------------------------------
sumx = sum(Xi);          % number of obs in each class
N = sum(sumx);          % total obs
repsum = sumx(ones(nfeatures, 1), :);

% priors_Pt, p(t) for each of t tasks
% assume Xi are task indicators, must be zero or one
priors_Pt = sumx ./ N;


%priors_Pt = sum(Xi) ./ sum(Xi(:));


% evidence, p(a) for each activation state S = {0, 1}

% Note: Only return this in 'full' output mode, because classification can be
% (is) made on the likelihood * p(t), which is proportional to p(t | a)
% also, this way, we can accumulate likelihood based on the dist. of all
% the data, and apply the priors_Pt only once.

% % % no need to save?
% % evidence_yes = sum(Y) ./ N;
% % evidence_no = 1 - evidence_yes;

% Likelihoods and posteriors: conditional probabilities
% -----------------------------------------------------

% likelihoods: p(activity | task)
xcount = Y' * Xi;       % activation count for each task
xnotcount = repsum - xcount;

% k is regularization param, biases towards 0.5 (2 classes) or 1/# classes
pa1_given_t = (xcount + k) ./ (nclasses*k + repsum); %(repmat(sumx, nfeatures, 1));

pa0_given_t = (xnotcount + k) ./ (nclasses*k + repsum); %(repmat(sumx, nfeatures, 1));


% % % Likelihood ratio
% % evidence_yes = (sum(Y) ./ nobs)';
% % evidence_yes = evidence_yes(:, ones(1, nclasses));
% % % add .1 to regularize...
% % lr = ( (pa1_given_t + .1) ./ (evidence_yes + .1) ) - 1;

%%% need to apply priors_Pt only once!! Not to each voxel. Don't need to
%%% return full posterior for each voxel.
% posteriors: p(task | activity)
if nargout > 3

    notY = ~Y;
    sumv = sum(Y);
    sumnotv = sum(notY);

    pt_given_act1 = xcount ./ (sumv(ones(2, 1), :)');
    pt_given_act0 = (notY' * Xi) ./ (sumnotv(ones(2, 1), :)');

    pt_given_act1(pt_given_act1 == 0) = .001;
    pt_given_act1(pt_given_act1 == 1) = .999;

    pt_given_act0(pt_given_act0 == 0) = .001;
    pt_given_act0(pt_given_act0 == 1) = .999;
    
    varargout{1} = pt_given_act1;
    varargout{2} = pt_given_act0;
end


if nargout > 5
    pa1_given_not_t = xnotcount ./ (N - repsum);
    
    pa1_given_not_t(pa1_given_not_t == 0) = .001;
    pa1_given_not_t(pa1_given_not_t == 1) = .999;
    
    varargout{3} = pa1_given_not_t;
end




% % % %     % priors_Pt, p(t) for each of t tasks
% % % %     % assume Xi are task indicators, must be zero or one
% % % %     priors_Pt = sum(Xi) ./ sum(sum(Xi));
% % % %
% % % %     % P(act | task), sum(active and task) / sum(task)
% % % %     % P(Y | Xi)
% % % %     % nfeatures x nclasses
% % % %     % assume act. is 0 or 1
% % % %     %fprintf('Calculating likelihood. ');
% % % %
% % % %     % We add 1 to implement "add 1" or Laplace smoothing
% % % %     % so that classes will never have estimates of exactly zero.
% % % %     % This avoids the problem of sparseness w/a limited training set.
% % % %     % and the line below is no longer needed
% % % %     const = 1;       % not sure if this is right.
% % % %     % but adding any constant means that the less frequent classes will
% % % %     % have higher estimates for pa_ if there are no activations; that means
% % % %     % that a study that activates where no other studies did will tend to
% % % %     % be classified as the less frequent task!  not advantageous...
% % % %     % p(A = yes | T = t)
% % % %     %pa_given_t = ((Xi' * Y)' + 1) ./ (repmat(sum(Xi), nfeatures, 1) + const);
% % % %     pa_given_t = ((Xi' * Y)') ./ (repmat(sum(Xi), nfeatures, 1));
% % % %
% % % %     % make sure no zeros, b/c we can't take log of 0; regularize
% % % %     % This is important because it tells us how to treat the case when we
% % % %     % have *no* activations of a particular class in a voxel
% % % %     % small values here mean that we think it's *very* unlikely.
% % % %
% % % %     % this is always < # obs for any class, but is equal across classes, so
% % % %     % pa_ will be equal for all classes if there is no evidence.
% % % %     pa_given_t(pa_given_t < eps) = 1 ./ nobs;
% % % %     pa_given_t(pa_given_t == 1) = 1 - (1 ./ nobs) ;     % we use this in test to get p(notA | T=t), so we need this line too
% % % %
% % % %     % scaling factor to get actual probabilities p(task | act)
% % % %     % this is the Evidence, the joint p(act1, act2, ... etc.)
% % % %     % as long as there are no features with no activations, this will not
% % % %     % blow up.
% % % %     % even if no zeros in overall dataset, may be 0's in xvalidation
% % % %
% % % %     % joint probability of act1, act2, act3 ... activation in voxel1, 2, 3, etc.
% % % %     %fprintf('Calculating evidence. ');
% % % %     pa = sum(Y) ./ nobs;
% % % %     pa(pa < eps) = eps;
% % % %
% % % %     log_evidence_act = log( pa ) ;
% % % %     log_evidence_notact = log ( 1 - pa );
% % % %     true_class = indic2condf(Xi);
% % % %
% % % %     bayes_model = struct('nobs', nobs, 'nfeatures', nfeatures, 'nclasses', nclasses, ...
% % % %         'whkeep', whkeep, 'whkeep_from_notempty', whkeep_from_notempty, ...
% % % %         'priors_Pt', priors_Pt, 'pa_given_t', pa_given_t, 'log_evidence_act', log_evidence_act, 'log_evidence_notact', log_evidence_notact, ...
% % % %         'Xi', Xi, 'true_class', true_class);

%fprintf('Done. \n');
