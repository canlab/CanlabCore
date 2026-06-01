classdef cv_splitter
    % cv_splitter  Cross-validation splitter (sklearn-style).
    %
    % Single class with a `type` discriminator, dispatched at split() time.
    % Static factory methods provide ergonomic constructors:
    %
    %   sp = cv_splitter.kfold(5)
    %   sp = cv_splitter.stratified_kfold(10)
    %   sp = cv_splitter.group_kfold(5)
    %   sp = cv_splitter.stratified_group_kfold(5)
    %   sp = cv_splitter.leave_one_group_out()
    %   sp = cv_splitter.repeated_kfold(5, 10)
    %   sp = cv_splitter.shuffle_split(10, 0.2)
    %   sp = cv_splitter.holdout(0.2)
    %   sp = cv_splitter.custom_partition(fold_ids)
    %
    % Use:
    %   splits = sp.split(X, Y, groups)         % all args
    %   splits = sp.split(X)                    % when not stratified / grouped
    %
    % Returns a 1 x nfolds struct array with logical .trIdx and .teIdx.
    %
    % All splitters set rng(obj.random_state) if random_state is set, for
    % reproducibility. The "groups" argument is the subject id vector for
    % leave-one-subject-out designs — its name and semantics mirror
    % sklearn's `groups`.
    %
    % :Example:
    % ::
    %   sp = cv_splitter.stratified_group_kfold(5);
    %   sp.random_state = 0;
    %   splits = sp.split(X, Y, subject_id);
    %   for f = 1:numel(splits)
    %       tr = splits(f).trIdx; te = splits(f).teIdx;
    %       % fit on X(tr,:), Y(tr); test on X(te,:), Y(te)
    %   end

    properties
        type             = '';   % 'kfold','stratified_kfold','group_kfold',
                                 % 'stratified_group_kfold','logo','holdout',
                                 % 'shuffle_split','repeated_kfold','custom'
        k                = 5;    % number of folds
        nrepeats         = 1;    % for repeated_kfold
        test_size        = 0.2;  % for holdout/shuffle_split (fraction in (0,1) or count >=1)
        n_splits         = 10;   % for shuffle_split
        random_state     = [];   % rng seed (empty = no reseed)
        custom_folds     = [];   % integer vector for custom_partition
    end

    methods
        % -----------------------------------------------------------------
        function obj = cv_splitter(type, varargin)
            if nargin == 0, return; end
            obj.type = type;
            p = inputParser; p.KeepUnmatched = false;
            addParameter(p, 'k',            obj.k);
            addParameter(p, 'nrepeats',     obj.nrepeats);
            addParameter(p, 'test_size',    obj.test_size);
            addParameter(p, 'n_splits',     obj.n_splits);
            addParameter(p, 'random_state', obj.random_state);
            addParameter(p, 'custom_folds', obj.custom_folds);
            parse(p, varargin{:});
            fn = fieldnames(p.Results);
            for i = 1:numel(fn), obj.(fn{i}) = p.Results.(fn{i}); end
        end

        % -----------------------------------------------------------------
        function splits = split(obj, X, Y, groups)
            if nargin < 4, groups = []; end
            if nargin < 3, Y      = []; end

            if ~isempty(obj.random_state)
                rng(obj.random_state);
            end

            if isnumeric(X) || islogical(X)
                n = size(X, 1);
            else
                error('cv_splitter:split:BadX', 'X must be a numeric or logical matrix.');
            end

            switch obj.type
                case 'kfold'
                    splits = cv_splitter.do_kfold(n, obj.k);
                case 'stratified_kfold'
                    if isempty(Y), error('stratified_kfold needs Y.'); end
                    splits = cv_splitter.do_stratified_kfold(Y, obj.k);
                case 'group_kfold'
                    if isempty(groups), error('group_kfold needs groups.'); end
                    splits = cv_splitter.do_group_kfold(groups, obj.k);
                case 'stratified_group_kfold'
                    if isempty(Y) || isempty(groups)
                        error('stratified_group_kfold needs Y and groups.');
                    end
                    splits = cv_splitter.do_stratified_group_kfold(Y, groups, obj.k);
                case 'logo'
                    if isempty(groups), error('leave_one_group_out needs groups.'); end
                    splits = cv_splitter.do_logo(groups);
                case 'holdout'
                    splits = cv_splitter.do_holdout(n, obj.test_size);
                case 'shuffle_split'
                    splits = cv_splitter.do_shuffle_split(n, obj.n_splits, obj.test_size);
                case 'repeated_kfold'
                    splits = cv_splitter.do_repeated_kfold(n, obj.k, obj.nrepeats);
                case 'custom'
                    if isempty(obj.custom_folds), error('custom needs custom_folds.'); end
                    splits = cv_splitter.do_custom(obj.custom_folds);
                otherwise
                    error('cv_splitter:UnknownType', ...
                        'Unknown splitter type: %s', obj.type);
            end
        end
    end


    methods (Static)
        % -------------------- Factory functions --------------------
        function s = kfold(k),                       s = cv_splitter('kfold', 'k', k); end
        function s = stratified_kfold(k),            s = cv_splitter('stratified_kfold', 'k', k); end
        function s = group_kfold(k),                 s = cv_splitter('group_kfold', 'k', k); end
        function s = stratified_group_kfold(k),      s = cv_splitter('stratified_group_kfold', 'k', k); end
        function s = leave_one_group_out(),          s = cv_splitter('logo'); end
        function s = repeated_kfold(k, nrepeats),    s = cv_splitter('repeated_kfold', 'k', k, 'nrepeats', nrepeats); end
        function s = shuffle_split(n_splits, test_size)
            s = cv_splitter('shuffle_split', 'n_splits', n_splits, 'test_size', test_size);
        end
        function s = holdout(test_size),             s = cv_splitter('holdout', 'test_size', test_size); end
        function s = custom_partition(fold_ids),     s = cv_splitter('custom', 'custom_folds', fold_ids); end

        % -------------------- Split implementations --------------------
        function splits = do_kfold(n, k)
            order = randperm(n)';
            fold_id = zeros(n, 1);
            fold_id(order) = mod((0:n-1)', k);
            splits = cv_splitter.fold_id_to_splits(fold_id);
        end

        function splits = do_stratified_kfold(Y, k)
            Y = Y(:); n = numel(Y);
            classes = unique(Y(~isnan(Y)));
            fold_id = zeros(n, 1);
            for c = 1:numel(classes)
                idx = find(Y == classes(c));
                idx = idx(randperm(numel(idx)));
                fold_id(idx) = mod((0:numel(idx)-1)', k);
            end
            splits = cv_splitter.fold_id_to_splits(fold_id);
        end

        function splits = do_group_kfold(groups, k)
            groups = groups(:);
            uniq = unique(groups);
            uniq = uniq(randperm(numel(uniq)));
            group_fold = mod((0:numel(uniq)-1)', k);
            fold_id = zeros(numel(groups), 1);
            for i = 1:numel(uniq)
                fold_id(groups == uniq(i)) = group_fold(i);
            end
            splits = cv_splitter.fold_id_to_splits(fold_id);
        end

        function splits = do_stratified_group_kfold(Y, groups, k)
            Y = Y(:); groups = groups(:);
            uniq = unique(groups);
            % Each group's "label" for stratification = mode of its Y values.
            % (Assumes within-group Y is mostly constant — typical for
            % subject-level designs.)
            group_label = zeros(numel(uniq), 1);
            for i = 1:numel(uniq)
                group_label(i) = mode(Y(groups == uniq(i)));
            end
            classes = unique(group_label);
            group_fold = zeros(numel(uniq), 1);
            for c = 1:numel(classes)
                idx = find(group_label == classes(c));
                idx = idx(randperm(numel(idx)));
                group_fold(idx) = mod((0:numel(idx)-1)', k);
            end
            fold_id = zeros(numel(groups), 1);
            for i = 1:numel(uniq)
                fold_id(groups == uniq(i)) = group_fold(i);
            end
            splits = cv_splitter.fold_id_to_splits(fold_id);
        end

        function splits = do_logo(groups)
            groups = groups(:);
            uniq = unique(groups);
            k = numel(uniq);
            splits = struct('trIdx', cell(1, k), 'teIdx', cell(1, k));
            for i = 1:k
                te = groups == uniq(i);
                splits(i).trIdx = ~te;
                splits(i).teIdx = te;
            end
        end

        function splits = do_holdout(n, test_size)
            n_test = cv_splitter.size_to_count(test_size, n);
            perm = randperm(n);
            te = false(n, 1); te(perm(1:n_test)) = true;
            splits.trIdx = ~te;
            splits.teIdx = te;
        end

        function splits = do_shuffle_split(n, n_splits, test_size)
            n_test = cv_splitter.size_to_count(test_size, n);
            splits = struct('trIdx', cell(1, n_splits), 'teIdx', cell(1, n_splits));
            for i = 1:n_splits
                perm = randperm(n);
                te = false(n, 1); te(perm(1:n_test)) = true;
                splits(i).trIdx = ~te;
                splits(i).teIdx = te;
            end
        end

        function splits = do_repeated_kfold(n, k, nrepeats)
            splits = struct('trIdx', cell(1, k*nrepeats), 'teIdx', cell(1, k*nrepeats));
            j = 0;
            for r = 1:nrepeats
                rep = cv_splitter.do_kfold(n, k);
                for f = 1:numel(rep)
                    j = j + 1;
                    splits(j) = rep(f);
                end
            end
        end

        function splits = do_custom(custom_folds)
            cf = custom_folds(:);
            uniq = unique(cf);
            k = numel(uniq);
            splits = struct('trIdx', cell(1, k), 'teIdx', cell(1, k));
            for i = 1:k
                te = cf == uniq(i);
                splits(i).trIdx = ~te;
                splits(i).teIdx = te;
            end
        end

        function splits = fold_id_to_splits(fold_id)
            fold_id = fold_id(:);
            uniq = unique(fold_id);
            k = numel(uniq);
            splits = struct('trIdx', cell(1, k), 'teIdx', cell(1, k));
            for i = 1:k
                te = fold_id == uniq(i);
                splits(i).trIdx = ~te;
                splits(i).teIdx = te;
            end
        end

        function c = size_to_count(test_size, n)
            % Convert a fraction-in-(0,1) or absolute count >=1 to a sample count.
            if test_size <= 0 || test_size >= n
                error('test_size %g not valid for n=%d', test_size, n);
            end
            if test_size < 1
                c = round(test_size * n);
            else
                c = round(test_size);
            end
            c = max(1, min(c, n - 1));
        end
    end

end
