function [b_star, gam_hat, Vgam_hat, gam_t, sigma2_b_est, subject_w] = RB_empirical_bayes_params(b, V)
%[b_star, gam_hat, Vgam_hat, gam_t, sigma2_b_est, subject_w] = RB_empirical_bayes_params(b, V)
%
% Stuff from Raudenbush & Bryk book
%
% Input: 
% b, betas from each subject
%    ** NOTE: ASSUMES b is regressors OF INTEREST
%    If nuisance covariates are added, V may be larger than b
%    and this program assumes the additional rows/cols of V represent
%    nuisance covariates and ignores them.
%
% V = xtxi * sigmasq;                           % Var/Cov mtx, Precision^-1, used in weighted est. and empirical bayes
%
% Fixed effects
% gam : pop. fixed effect
% gam_hat, est of gam
% W : 2nd-level design matrix for gam, would contain indicators for group
% membership
%     is 1st level pred x groups, F in 3.26 (p. 43) is groups; gam is F x
%     1, one param for each group ?? could be diff W for each subj
% u : random coefficients, distributed N(0, T),
% T is residual p x p var/cov matrix for 2nd-level effects
% matrix for 2nd level; eq. 3.26, p 43.
% D : delta, dispersion matrix, inverse of precision matrix
%     D = T + V{i}
% Precis, precision matrix (inv of total b + w cov matrix), D^-1, 3.29
% Z : Columns of X with random effects
% Vy : "marginal" Var/Cov matrix of y
 
Nsubj = length(V);

% Get dimensionality of predictor cov matrix
i = 1; Npred = 0; 
while Npred == 0, Npred = size(V{i}, 2); end
 
% Handle empty matrices (if no cov data for some bs)
% also handle case where V has more cols/rows than b has rows (predictors)
nregs = size(b, 1);

for i = 1:Nsubj
    if isempty(V{i}), V{i} = Inf .* eye(Npred); 

    elseif size(V{i}, 1) > nregs, V{i} = V{i}(1:nregs, 1:nregs); 
        
    end
end
    
    
T = nancov(b');  % between-subjects cov matrix; TOR's naive estimate (not R & B); could at least subtract expected w/i ss contrib?
for i = 1:Nsubj, Dind{i} = T + V{i}; Precis{i} = inv(Dind{i});  end
 
W = eye(nregs);  % 2nd-level design matrix; intercept only = simple
 

% Estimate fixed effects (gammas), Eq. 3.31 of R & B
% ------------------------------------------------------------
% This is WLS estimate, whereas naive est. would be mean(b, 2)
% eq 3.31, use blk diag, simpler...
% Note: if W = I, inv(blkdiag(Precis{:})) == blkdiag(Dind{:})
G1 = zeros(size(W)); for i = 1:Nsubj, G1 = G1 + Precis{i}; end
G2 = zeros(nregs, 1); for i = 1:Nsubj, G2 = G2 + Precis{i} * b(:, i); end
gam_hat = inv(G1) * G2;
% ------------------------------------------------------------

% Inference on fixed effects (gammas)
% ------------------------------------------------------------
% inv(G1) is Var(gam-hat)
Vgam_hat = inv(G1);
gam_t = gam_hat ./ (diag(Vgam_hat) .^ .5);
% ------------------------------------------------------------


% Estimate individual slopes (betas, fixed + random)
% ------------------------------------------------------------
% Now we have 3 choices for betas: the individual b estimates, the fixed
% effect estimates, or an Empirical Bayes combination of the two
% L is "optimal" weighting matrix for each subject
% R & B eq. 3.56
I = eye(size(Precis{1}));
for i = 1:Nsubj
    L{i} = T * Precis{i};   % T(T + V{i})^-1
    b_star(:, i) = L{i} * b(:, i) + (I - L{i}) * W * gam_hat;  % weighted linear combo  of indiv and fixed fx est.
end
 
% Tor made up the stuff below...develop...
% ---------------------------------------------
Vb_star = cov(b_star');
 
if Nsubj == 1 
    % can't calculate; return NaN
    sigma2_b_est = NaN * ones(Npred, 1);
else
    sigma2_b_est = diag(Vb_star);  % this is shrunk from original est.
end

for i = 1:Nsubj, sigma2_tot_est(:, i) = sigma2_b_est + diag(V{i}); end
subject_w = (1 ./  sigma2_tot_est)';
 
subject_w = subject_w ./ repmat(sum(subject_w), Nsubj, 1);
 
 
end
