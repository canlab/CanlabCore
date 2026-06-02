function obj = permutation_test(obj, X, Y, varargin)
% permutation_test  Build a null distribution by shuffling Y and re-running CV.
%
% :Usage:
% ::
%     pm = permutation_test(pm, X, Y);
%     pm = permutation_test(pm, X, Y, 'nperm', 5000);
%     pm = permutation_test(pm, X, Y, 'groups', subject_id);
%
% Pipeline:
%   1. Run crossval(obj, X, Y) once to get the observed score.
%   2. For each of nperm permutations:
%        - shuffle Y (preserving groups if 'groups' provided —
%          subjects are permuted as units)
%        - re-run crossval on (X, Y_perm) with the same cv splitter
%          and scorer
%        - record the null cv score
%   3. p-value = (1 + #null >= observed) / (1 + nperm), with
%      direction flipped for less-is-better scorers.
%
% After:
%   obj.permutation_results.observed     observed cv score
%   obj.permutation_results.null_scores  [nperm x 1]
%   obj.permutation_results.p_value      one-sided permutation p
%   obj.permutation_results.scorer       scorer name used
%   obj.permutation_results.nperm

    pi = inputParser; pi.KeepUnmatched = true;
    addParameter(pi, 'nperm',   []);
    addParameter(pi, 'groups',  []);
    addParameter(pi, 'verbose', true);
    parse(pi, varargin{:});

    nperm = pi.Results.nperm;
    if isempty(nperm)
        nperm = obj.nperm;
        if isempty(nperm), nperm = 1000; end
    end
    groups  = pi.Results.groups;
    verbose = pi.Results.verbose;

    if ~isempty(obj.random_state), rng(obj.random_state); end

    % 1. Observed score via crossval.
    obs_obj   = crossval(obj, X, Y, 'groups', groups);
    scorer    = obs_obj.scorer;
    metric    = scorer.name;
    obs_score = obs_obj.error_metrics.(metric).value;

    % 2. Null distribution.
    null_scores = nan(nperm, 1);
    n = numel(Y); Y = Y(:);
    if verbose, fprintf('permutation_test: %d permutations', nperm); end
    every = max(1, ceil(nperm / 20));

    for i = 1:nperm
        if verbose && mod(i, every) == 0, fprintf('.'); end

        if isempty(groups)
            % No grouping: free observation-level permutation.
            Y_perm = Y(randperm(n));
        else
            uniq = unique(groups);
            % Detect whether Y is constant within each group. If yes,
            % the design is "between-subjects" (each subject has one
            % class) and we should shuffle Y at the SUBJECT level. If
            % no, the design is "within-subjects" (each subject has
            % multiple classes; e.g. DPSP Hot+Warm per subject) and we
            % shuffle Y WITHIN each subject — the correct permutation
            % for paired data.
            y_within_const = arrayfun(@(g) ...
                numel(unique(Y(groups == g))) == 1, uniq);
            if all(y_within_const)
                % Subject-level shuffle.
                group_Y = arrayfun(@(g) Y(find(groups == g, 1, 'first')), uniq);
                new_group_Y = group_Y(randperm(numel(uniq)));
                Y_perm = Y;
                for g = 1:numel(uniq)
                    Y_perm(groups == uniq(g)) = new_group_Y(g);
                end
            else
                % Within-subject shuffle: permute Y for each subject's
                % observations independently.
                Y_perm = Y;
                for g = 1:numel(uniq)
                    ix = find(groups == uniq(g));
                    Y_perm(ix) = Y(ix(randperm(numel(ix))));
                end
            end
        end

        try
            null_obj = crossval(clone(obj), X, Y_perm, 'groups', groups);
            null_scores(i) = null_obj.error_metrics.(metric).value;
        catch
            % skip
        end
    end
    if verbose, fprintf(' done\n'); end

    % 3. p-value
    valid = ~isnan(null_scores);
    if scorer.greater_is_better
        n_extreme = sum(null_scores(valid) >= obs_score);
    else
        n_extreme = sum(null_scores(valid) <= obs_score);
    end
    pval = (1 + n_extreme) / (1 + sum(valid));

    obj = obs_obj;  % start from the populated obs object
    obj.permutation_results.observed    = obs_score;
    obj.permutation_results.null_scores = null_scores;
    obj.permutation_results.p_value     = pval;
    obj.permutation_results.scorer      = metric;
    obj.permutation_results.nperm       = nperm;

    obj.history{end+1, 1} = sprintf( ...
        'permutation_test(%s): observed=%.3f, p=%.4f (nperm=%d)', ...
        metric, obs_score, pval, nperm);
end
