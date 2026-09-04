function R = hrf_mediation_analyze(X, M, Y, varargin)
%HRF_MEDIATION_ANALYZE Trial-level (multilevel) mediation X -> M -> Y.
%
% Tests whether a mediator M carries the effect of X on Y at the trial level:
%   a  : X -> M
%   b  : M -> Y  (controlling for X)
%   c  : X -> Y  (total)
%   c' : X -> Y  (direct, controlling for M)
%   ab : a*b     (INDIRECT / mediated effect)   -- the quantity of interest
%
% This is the IO-free compute core of the HRF mediation method (the driver
% hrf_causality_mediation builds X/M/Y from deconvolved node amplitudes and/or
% the BIDS events, e.g. stimulus -> region -> rating). It is what makes
% "does region/signature M mediate stimulus -> behaviour" answerable on top of
% the same deconvolution used for the Granger analysis.
%
% Multilevel (>=2 subjects): the five paths are estimated per subject by OLS,
% then summarized across subjects with a one-sample t-test (the standard
% two-level mediation summary). Single subject / pooled: the indirect effect
% ab is tested by a nonparametric bootstrap over trials.
%
% :Usage:
% ::
%     R = hrf_mediation_analyze(Xc, Mc, Yc, 'Names', {'stim','NPS','rating'})
%     R = hrf_mediation_analyze(x, m, y, 'Nboot', 10000)        % single-level
%
% :Inputs:
%   **X, M, Y:** trial-level data. Either a cell{1..nSubj} of column vectors
%             (one per subject; same length within a subject) for multilevel
%             mediation, or a single numeric vector for one subject / pooled.
%
% :Optional Inputs:
%   **'Names':**       {Xname, Mname, Yname} for labelling. Default generic.
%   **'Covariates':**  cell{1..nSubj} of [nTrial x q] nuisance regressors
%                      partialled from all paths (or a single matrix). Default none.
%   **'Standardize':** z-score X, M, Y within subject first (default true), so
%                      paths are in comparable (standardized) units.
%   **'Nboot':**       bootstrap samples for the single-level indirect test.
%                      Default 5000.
%   **'Verbose'/'doverbose':** print a summary (default true).
%
% :Output:
%   **R:** struct with one sub-struct per path (a, b, cp, c, ab), each holding
%          .est (group/point estimate), .se, .t, .p, and .persubj (the
%          per-subject estimates, multilevel). Plus .prop_mediated (ab/c),
%          .n_subj, .n_trials, .names, .level ('multilevel'|'single').
%
% :Examples:
% ::
%     % synthetic: X->M->Y with a=0.6, b=0.5
%     ns=12; X=cell(1,ns); M=X; Y=X;
%     for s=1:ns, x=randn(40,1); m=0.6*x+randn(40,1); y=0.5*m+0.2*x+randn(40,1);
%         X{s}=x; M{s}=m; Y{s}=y; end
%     R = hrf_mediation_analyze(X,M,Y,'Names',{'stim','region','rating'});
%     R.ab.est, R.ab.p          % indirect effect ~0.3, significant
%
% See also: hrf_causality_mediation, hrf_causality, mediation (MediationToolbox).

p = inputParser;
p.addRequired('X');
p.addRequired('M');
p.addRequired('Y');
p.addParameter('Names', {}, @(x) iscell(x) || isstring(x));
p.addParameter('Covariates', {}, @(x) iscell(x) || isnumeric(x));
p.addParameter('Standardize', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Nboot', 5000, @(x) isscalar(x) && x >= 100);
% path significance at the group (multilevel) level: 'none' (default,
% parametric one-sample t) or 'permutation' (per-path sign-flip; exact at
% n<=13, robust at small n).
p.addParameter('Correction', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('doverbose', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.parse(X, M, Y, varargin{:});
opts = p.Results;
verbose = logical(opts.Verbose);
if ~isempty(opts.doverbose), verbose = logical(opts.doverbose); end

Xc = local_to_cell(X); Mc = local_to_cell(M); Yc = local_to_cell(Y);
nSubj = numel(Xc);
if numel(Mc) ~= nSubj || numel(Yc) ~= nSubj
    error('hrf_mediation_analyze:Mismatch', 'X, M, Y must have the same number of subjects.');
end
Cov = local_cov_cell(opts.Covariates, nSubj);
names = local_names(opts.Names);

% per-subject paths
paths = nan(nSubj, 5);   % [a b cp c ab]
ntrials = zeros(nSubj, 1);
for s = 1:nSubj
    [x, m, y, cv, ok] = local_prep(Xc{s}, Mc{s}, Yc{s}, Cov{s}, logical(opts.Standardize));
    if ~ok, continue; end
    ntrials(s) = numel(x);
    paths(s, :) = local_paths(x, m, y, cv);
end
valid = all(isfinite(paths), 2);
paths = paths(valid, :); nValid = size(paths, 1);
if nValid == 0, error('hrf_mediation_analyze:NoData', 'No subject had usable trials.'); end

if nValid >= 2
    level = 'multilevel';
    R = local_group_struct(paths, names, nValid, opts.Correction, opts.Nboot);
else
    level = 'single';
    [x, m, y, cv] = local_prep(Xc{1}, Mc{1}, Yc{1}, Cov{1}, logical(opts.Standardize));
    R = local_boot_struct(x, m, y, cv, opts.Nboot, names);
end
R.level = level;
R.n_subj = nValid;
R.n_trials = sum(ntrials);
R.names = names;

if verbose
    fprintf('hrf_mediation_analyze [%s]: %s -> %s -> %s | n=%d %s, %d trials\n', ...
        level, names{1}, names{2}, names{3}, nValid, ...
        local_tern(strcmp(level, 'multilevel'), 'subjects', 'subject'), R.n_trials);
    fprintf('  a=%+.3f(p=%.2g)  b=%+.3f(p=%.2g)  c=%+.3f  c''=%+.3f  | INDIRECT ab=%+.3f (p=%.2g)\n', ...
        R.a.est, R.a.p, R.b.est, R.b.p, R.c.est, R.cp.est, R.ab.est, R.ab.p);
    if isfinite(R.prop_mediated)
        fprintf('  proportion mediated = %.0f%%\n', 100 * R.prop_mediated);
    end
end
end


% =========================================================================
function ab_paths = local_paths(x, m, y, cov)
% OLS path coefficients for one subject. cov is [n x q] (may be empty).
n = numel(x);
Da = [ones(n, 1), x, cov];          ba = Da \ m; a = ba(2);
Db = [ones(n, 1), x, m, cov];       bb = Db \ y; cp = bb(2); b = bb(3);
Dc = [ones(n, 1), x, cov];          bc = Dc \ y; c = bc(2);
ab_paths = [a, b, cp, c, a * b];
end


function R = local_group_struct(paths, names, n, correction, nperm)
% Group stats for each path across subjects. Default parametric one-sample t;
% Correction='permutation' uses a sign-flip test PER PATH (exact at n<=13,
% robust at small n) -- no cross-path correction (the 5 paths are distinct
% questions, not one family).
lab = {'a', 'b', 'cp', 'c', 'ab'};
useperm = any(strcmpi(char(correction), {'permutation', 'perm', 'maxt'}));
R = struct();
for k = 1:5
    v = paths(:, k);
    mu = mean(v, 'omitnan');
    se = std(v, 0, 'omitnan') / sqrt(n);
    t = mu / se; if se == 0, t = 0; end
    if useperm
        G = hrf_group_stats(reshape(v, 1, n), 'Correction', 'permutation', 'Nperm', nperm);
        pv = G.p_corr;                    % single cell -> the sign-flip p
    else
        pv = 2 * (1 - local_tcdf(abs(t), n - 1));
    end
    R.(lab{k}) = struct('est', mu, 'se', se, 't', t, 'p', pv, 'persubj', v);
end
R.prop_mediated = R.ab.est / R.c.est;
if ~isfinite(R.prop_mediated), R.prop_mediated = NaN; end
R.names = names;
end


function R = local_boot_struct(x, m, y, cov, nboot, names)
% Single-level: point paths + bootstrap CI/p for the indirect effect.
pt = local_paths(x, m, y, cov);
lab = {'a', 'b', 'cp', 'c', 'ab'};
n = numel(x);
boot = nan(nboot, 5);
for bms = 1:nboot
    idx = local_bootidx(n, bms);          % deterministic per-iter resample
    cb = cov; if ~isempty(cb), cb = cb(idx, :); end
    boot(bms, :) = local_paths(x(idx), m(idx), y(idx), cb);
end
R = struct();
for k = 1:5
    bk = boot(:, k);
    se = std(bk, 'omitnan');
    frac_pos = mean(bk > 0); pv = 2 * min(frac_pos, 1 - frac_pos); pv = max(pv, 1 / nboot);
    R.(lab{k}) = struct('est', pt(k), 'se', se, 't', pt(k) / se, 'p', pv, ...
        'ci', quantile(bk, [0.025 0.975]), 'persubj', []);
end
R.prop_mediated = pt(5) / pt(4);
if ~isfinite(R.prop_mediated), R.prop_mediated = NaN; end
R.names = names;
end


% =========================================================================
function [x, m, y, cv, ok] = local_prep(x, m, y, cv, standardize)
x = double(x(:)); m = double(m(:)); y = double(y(:));
n = min([numel(x), numel(m), numel(y)]);
x = x(1:n); m = m(1:n); y = y(1:n);
if ~isempty(cv), cv = double(cv(1:min(size(cv, 1), n), :)); end
good = isfinite(x) & isfinite(m) & isfinite(y);
if ~isempty(cv), good = good & all(isfinite(cv), 2); end
x = x(good); m = m(good); y = y(good);
if ~isempty(cv), cv = cv(good, :); end
ok = numel(x) >= 5 && std(x) > 0 && std(m) > 0 && std(y) > 0;
if ok && standardize
    x = local_z(x); m = local_z(m); y = local_z(y);
    if ~isempty(cv), cv = local_z(cv); end
end
end

function z = local_z(v)
mu = mean(v, 1); sd = std(v, 0, 1); sd(sd == 0) = 1;
z = (v - mu) ./ sd;
end

function c = local_to_cell(v)
if iscell(v), c = v(:)'; else, c = {v}; end
end

function C = local_cov_cell(cov, n)
if isempty(cov), C = repmat({[]}, 1, n); return; end
if iscell(cov), C = cov(:)'; if numel(C) ~= n, C = repmat(C(1), 1, n); end
else, C = repmat({cov}, 1, n); end
end

function names = local_names(nm)
def = {'X', 'M', 'Y'};
if isempty(nm), names = def; return; end
nm = cellstr(string(nm));
names = def; names(1:min(3, numel(nm))) = nm(1:min(3, numel(nm)));
end

function idx = local_bootidx(n, seed)
% Deterministic resample with replacement (avoids rand for reproducibility).
a = mod(1103515245 * seed + 12345, 2^31);
idx = zeros(n, 1);
for i = 1:n
    a = mod(1103515245 * a + 12345, 2^31);
    idx(i) = 1 + floor((a / 2^31) * n);
end
end

function pcdf = local_tcdf(tval, df)
x = df ./ (df + tval .^ 2);
pcdf = 1 - 0.5 * betainc(x, df / 2, 0.5);
end

function s = local_tern(c, a, b)
if c, s = a; else, s = b; end
end
