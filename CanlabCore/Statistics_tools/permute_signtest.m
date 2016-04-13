function [p, Z, xbar, permsign, pmean] = permute_signtest(data, nperms, w, permsign)
% One-sample t-test against zero on each column of data
% Using weighted sign permutation test
%
% :Usage:
% ::
%
%     [p, Z, xbar, permsign, pmean] = permute_signtest(data, nperms, [w], [permsign])
%
% p-values are 2-tailed
%
% :Inputs:
%
%   **data:**
%        n x k matrix of k vectors to test against zero
%
%   **nperms:**
%        number of permutations
%
%   **w:**
%        n x k matrix of weights
%
%        weights for each column should sum to 1
%
%        default (for empty or missing input) is equal weights, i.e., 1/n
%
%   **permindx:**
%        nperms x n matrix of exact permutation indices
%        default (for empty or missing input) generates nperms permutations
%
% :Example: Generate 5-vector dataset and test, first generating perms, then
%          using already-generated ones
% ::
%
%    tic, [p, z, xbar, permsign] = permute_signtest(x, 1000); toc
%    tic, [p2, z2, xbar, permsign] = permute_signtest(x, [], [], permsign); toc
%
%
% Example: Generate fake data and simulate false positive rate
% ::
%
%    permsign = [];  % initialize to empty for first iteration
%    nperms = 2000; nreps = 5000;
%    tic, for i = 1:nreps
%    x = randn(20,1); stat = mean(x);
%    [p(i), z(i), xbar, permsign] = permute_signtest(x, nperms, [], permsign);
%    if mod(i,10) == 0,fprintf(1,'%03d ',i); end
%    end
%    fprintf(1,'\n');
%
% Example: Generate fake data and simulate false positive rate
% This uses the column-wise capabilities of this function and is much
% faster
% % ::
%
%    nperms = 2000; nreps = 5000; nsubj = 15;
%    x = randn(nsubj, nreps);
%    [p, z, xbar, permsign] = permute_signtest(x, nperms);
%
% ..
%    Tor Wager, march 2007
% ..


[n, k] = size(data); % observations x columns

if nargin < 4 || isempty(permsign)
    [permsign, nperms] = setup_signperms(n, nperms); % rows are permutations
else
    nperms = size(permsign, 2);
end

% this would work for permuting rows, i.e., for correlation test
% % if nargin < 4 || isempty(permindx)
% %     permindx = permute_setupperms(n, nperms); % rows are permutations
% % else
% %     nperms = size(permindx, 1);
% % end

if nargin < 3 || isempty(w)
    w = ones(n, k) ./ n;
end

pmean = zeros(nperms, k);

wmean = @(data, w) diag(w'*data)';

% correct permutation

xbar = wmean(data, w);


% sign permutations
for i = 1:nperms

    signs = repmat(permsign(:,i), 1, k);

    pdata = data .* signs;

    % permutation-distribution means
    pmean(i,:)  = wmean(pdata, w);

end


xbarmtx = repmat(xbar, nperms, 1);
p = min( [sum(pmean <= xbarmtx); sum(pmean >= xbarmtx)] ) ./ nperms;

p = 2 .* p; % two-tailed

% make sure p is not 1 or 0 due to limited bootstrap samples
p = max(p, 1 ./ nperms);
p = min(p, 1 - 1./nperms);

% adjust signs of z-scores to reflect direction of effect
Z = abs(norminv(p ./ 2)) .* sign(xbar);

end



function [permsign, nperms] = setup_signperms(n, nperms)
% get all possible combinations of [1 -1] signs for n subjects
% 2^n possible combos

% re-set rand number generator, which is used by randperm
rand('twister',sum(100*clock))

% %   t1 = clock;
% %   fprintf(1, 'Setting up perms ');

% can only have as many perms as unique values.
% if we have *all* unique perms, this is an exact test
nperms = min(nperms, 2^n);

% columns are all possible combos

if 2^n <= 32768
    % feasible to get *all* sign permutations; n = 15 or less

    vals = cell(1,n);
    [vals{:}] = deal([1 -1]);

    permsign = combvec(vals{:});

    % omit correct permutation
    permsign(:,1) = [];

    nc = size(permsign,2);

    if isempty(nperms)
        % do nothing besides set nperms
        nperms = nc;

    elseif nperms < nc
        % pick random subset
        p = randperm(nc);
        permsign = permsign(:, p(1:nperms));

    elseif nperms > nc
        % % %       fprintf(1,'Warning: only %3.0f perms available. Min p-value is %3.6f\n', nc, 1./nc)
        nperms = nc;

    end

else
    % too many perms of n [-1 1] signs
    % generate random subset
    n_to_start_with = min(2 * nperms, 2 * 2^n);   % generate a finite number greater than needed; we want unique perms only

    permsign = -1 + 2 * round(rand(n, n_to_start_with));

    permsign = unique(permsign', 'rows')';

    % Fill in add'l perms, if any needed
    % allow some repeats; less precise, but will help in generating perms fast
    n_needed_now = nperms - size(permsign, 2);
    permsign = [permsign  -1 + 2 * round(rand(n, n_needed_now))];
end

% % %   fprintf(1, 'Done in %3.0f s', etime(clock, t1));

end




