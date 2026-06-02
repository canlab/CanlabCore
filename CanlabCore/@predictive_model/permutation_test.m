function obj = permutation_test(obj, X, Y, varargin)
% permutation_test  Build a null distribution by shuffling Y and re-running CV.
%
% :Usage:
% ::
%     pm = permutation_test(pm, X, Y);
%     pm = permutation_test(pm, X, Y, 'nperm', 5000);
%     pm = permutation_test(pm, X, Y, 'groups', subject_id);
%     pm = permutation_test(pm, X, Y, 'groups', id, 'permutation', 'within_subjects');
%
% Pipeline:
%   1. Run crossval(obj, X, Y) once to get the observed score.
%   2. For each of nperm permutations: shuffle Y per the chosen
%      permutation scheme, re-run crossval, record the null cv score.
%   3. p-value = (1 + #null at-least-as-extreme) / (1 + nperm), with
%      direction flipped for less-is-better scorers.
%
% Permutation schemes ('permutation' name/value):
%
%   'auto'             (default) Detect from the data:
%                      * no groups               -> 'free'
%                      * groups + Y constant per group -> 'between_subjects'
%                      * groups + Y varies per group   -> 'within_subjects'
%
%   'free'             Observation-level shuffle of Y. Valid only if
%                      observations are truly independent. For grouped
%                      data this breaks within-subject correlation and
%                      generally INFLATES false positives — use the
%                      group-aware schemes instead.
%
%   'between_subjects' Reassign each subject's Y (uniform across that
%                      subject's observations) to another subject's Y,
%                      drawn at random without replacement. The
%                      appropriate null for between-subjects designs
%                      (each subject contributes one class). Errors
%                      with a warning if Y varies WITHIN a subject —
%                      that means the design is paired and this scheme
%                      can't represent the right null.
%
%   'within_subjects'  Permute Y INDEPENDENTLY within each subject's
%                      observations. The GOLD STANDARD null for paired
%                      / within-subject designs (e.g. each subject
%                      contributes both classes): preserves the
%                      subject-level pattern (no subject confound),
%                      breaks only the class-vs-brain mapping that the
%                      classifier is supposed to be exploiting. Warns
%                      if Y is constant within every subject — the
%                      permutation would be a no-op.
%
% After:
%   obj.permutation_results.observed              observed cv score
%   obj.permutation_results.null_scores           [nperm x 1]
%   obj.permutation_results.p_value               permutation p
%   obj.permutation_results.scorer                scorer name used
%   obj.permutation_results.nperm
%   obj.permutation_results.permutation           scheme actually used
%                                                 (auto resolves to one
%                                                 of the three concrete
%                                                 names below)
%   obj.permutation_results.permutation_descrip   one-line description

    pi = inputParser; pi.KeepUnmatched = true;
    addParameter(pi, 'nperm',       []);
    addParameter(pi, 'groups',      []);
    addParameter(pi, 'permutation', 'auto');
    addParameter(pi, 'verbose',     true);
    parse(pi, varargin{:});

    nperm = pi.Results.nperm;
    if isempty(nperm)
        nperm = obj.nperm;
        if isempty(nperm), nperm = 1000; end
    end
    groups       = pi.Results.groups;
    perm_request = lower(char(pi.Results.permutation));
    verbose      = pi.Results.verbose;

    if ~isempty(obj.random_state), rng(obj.random_state); end

    Y = Y(:);
    n = numel(Y);

    % Resolve the permutation scheme once, before the loop. Validate
    % forced choices and warn on degenerate combinations.
    if isempty(groups)
        % Without groups, only 'free' is meaningful.
        if ~ismember(perm_request, {'auto', 'free'})
            warning('predictive_model:permutation_test:NoGroups', ...
                ['You requested permutation=''%s'' but did not pass ''groups''. ' ...
                 'Falling back to ''free'' observation-level shuffle.'], perm_request);
        end
        perm_scheme  = 'free';
        perm_descrip = 'free observation-level shuffle (no groups provided)';
        uniq         = [];
        y_const      = [];
    else
        uniq    = unique(groups);
        y_const = arrayfun(@(g) numel(unique(Y(groups == g))) == 1, uniq);
        all_const  = all(y_const);
        any_varies = any(~y_const);

        switch perm_request
            case 'auto'
                if all_const
                    perm_scheme  = 'between_subjects';
                    perm_descrip = 'auto -> between_subjects (Y constant within each group)';
                else
                    perm_scheme  = 'within_subjects';
                    perm_descrip = 'auto -> within_subjects (Y varies within at least one group)';
                end

            case 'free'
                perm_scheme  = 'free';
                perm_descrip = 'free observation-level shuffle (groups ignored; forced by user)';

            case {'between_subjects', 'between'}
                perm_scheme = 'between_subjects';
                if any_varies
                    warning('predictive_model:permutation_test:InvalidForce', ...
                        ['permutation=''between_subjects'' is invalid: Y varies WITHIN ' ...
                         'at least one group (paired design). This scheme assumes Y is constant ' ...
                         'within each group. Within-subject variation will be lost when Y is ' ...
                         'broadcast back from a single per-group value (we use the first ' ...
                         'observation''s Y per group). Did you mean ''within_subjects''?']);
                end
                perm_descrip = 'between_subjects (forced by user)';

            case {'within_subjects', 'within', 'paired'}
                perm_scheme = 'within_subjects';
                if all_const
                    warning('predictive_model:permutation_test:DegenerateForce', ...
                        ['permutation=''within_subjects'' is degenerate: Y is constant within ' ...
                         'every group. The within-subject shuffle is a no-op and the null ' ...
                         'distribution will be identical to the observed score. Did you mean ' ...
                         '''between_subjects''?']);
                end
                perm_descrip = 'within_subjects (paired-design gold standard; forced by user)';

            otherwise
                error('predictive_model:permutation_test:UnknownScheme', ...
                    ['Unknown permutation=''%s''. Valid options: ''auto'', ''free'', ' ...
                     '''between_subjects'', ''within_subjects''.'], perm_request);
        end
    end

    % 1. Observed score via crossval.
    obs_obj   = crossval(obj, X, Y, 'groups', groups);
    scorer    = obs_obj.scorer;
    metric    = scorer.name;
    obs_score = obs_obj.error_metrics.(metric).value;

    % 2. Null distribution.
    null_scores = nan(nperm, 1);
    if verbose
        fprintf('permutation_test: %d perms (%s)', nperm, perm_scheme);
    end
    every = max(1, ceil(nperm / 20));

    for i = 1:nperm
        if verbose && mod(i, every) == 0, fprintf('.'); end

        switch perm_scheme
            case 'free'
                Y_perm = Y(randperm(n));

            case 'between_subjects'
                % One Y per group; reassign across groups at random.
                group_Y     = arrayfun(@(g) Y(find(groups == g, 1, 'first')), uniq);
                new_group_Y = group_Y(randperm(numel(uniq)));
                Y_perm = Y;
                for g = 1:numel(uniq)
                    Y_perm(groups == uniq(g)) = new_group_Y(g);
                end

            case 'within_subjects'
                Y_perm = Y;
                for g = 1:numel(uniq)
                    ix = find(groups == uniq(g));
                    Y_perm(ix) = Y(ix(randperm(numel(ix))));
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
    obj.permutation_results.observed             = obs_score;
    obj.permutation_results.null_scores          = null_scores;
    obj.permutation_results.p_value              = pval;
    obj.permutation_results.scorer               = metric;
    obj.permutation_results.nperm                = nperm;
    obj.permutation_results.permutation          = perm_scheme;
    obj.permutation_results.permutation_descrip  = perm_descrip;

    obj.history{end+1, 1} = sprintf( ...
        'permutation_test(%s, %s): observed=%.3f, p=%.4f (nperm=%d)', ...
        metric, perm_scheme, obs_score, pval, nperm);
end
