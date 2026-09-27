% rsm: Container for Representational Similarity / Dissimilarity Matrices.
%
% -------------------------------------------------------------------------
% Features and philosophy
% -------------------------------------------------------------------------
%
% @rsm is a value-class container for similarity (RSM) or dissimilarity (RDM)
% matrices with attached condition metadata, optional replicate axis
% (subjects, sessions, runs, folds), and provenance.
%
% It is the canonical output of fmri_data.compute_rsm and the canonical input
% to RSA inference methods (cells, contrasts, ttest_contrasts, reliability,
% drift, rsa_lm, rsa_lme, compare_models, compare, rsa_parcelwise,
% searchlight_rsa).
%
% Value semantics: every method returns a new rsm rather than mutating in
% place. R2 = R.fisher_z() does not modify R.
%
% Soft dependency on the Kriegeskorte rsatoolbox: methods probe via
% @rsm/private/probe_rsatoolbox.m and fall back to stock implementations
% when rsa.* is not on the path.
%
% -------------------------------------------------------------------------
% Properties
% -------------------------------------------------------------------------
%
% dat               [k x k] or [k x k x N] numeric matrix
%                   Similarity (RSM) or dissimilarity (RDM) values.
%                   3rd dim N indexes replicates whose meaning is given by
%                   the .level and .replicate_table properties.
%
% is_dissimilarity  logical
%                   true  -> .dat is an RDM (lower = more similar)
%                   false -> .dat is an RSM (higher = more similar)
%
% metric            char/string
%                   Name of the distance/similarity metric used to construct
%                   .dat. One of {'correlation','spearman','cosine','euclidean',
%                   'seuclidean','mahalanobis','crossnobis','design',
%                   'categorical','distance','none'} or any user-supplied string.
%
% labels            {k x 1} cellstr
%                   Condition / row-column labels. Length must equal size(dat,1).
%
% metadata_table    table with k rows (or empty)
%                   Per-condition attributes. One row per label. Optional.
%
% groupings         struct
%                   Map name -> indices into 1:k. Used by rsm.cells, contrast,
%                   plot 'block_borders_by', etc. Example:
%                       groupings.hot     = 1:8;
%                       groupings.warm    = 9:16;
%                       groupings.imagine = 17:24;
%
% level             char/string
%                   What each slice along dim 3 represents. One of:
%                   {'subject','session','run','collapsed','group','image',
%                    'model_stack','none'}.
%
% replicate_table   table with N rows (or empty)
%                   One row per slice along dim 3. Self-describes the
%                   replicate axis. Typical columns: subject_id,
%                   session_number, run_number, fold.
%
% whitened          struct
%                   Provenance of any whitening applied.
%                       .level        'none'|'within_subject'|'across_subject'|'session_difference'
%                       .method       'none'|'covdiag'|'diag'|'custom'
%                       .shrinkage    [] or scalar in [0,1]
%
% source            char/string
%                   Where this rsm came from. e.g., 'fmri_data',
%                   'parcel:<name>', 'searchlight:<vox>', 'custom', 'design'.
%
% history           cellstr
%                   Chronological log of operations applied to this rsm.
%
% additional_info   struct
%                   Free-form bag for paradigm-specific metadata.
%
% -------------------------------------------------------------------------
% Construction
% -------------------------------------------------------------------------
%
% R = rsm();                                          % empty
% R = rsm(dat);                                       % bare matrix, defaults filled
% R = rsm(dat, 'labels', {'A','B','C'}, ...
%              'metric', 'correlation', ...
%              'is_dissimilarity', false);
%
% R = rsm.from_categorical(cat_vec)                   % static — same-vs-different RDM
% R = rsm.from_metadata_distance(table, col, ...)     % static — continuous-distance RDM
% R = rsm.from_design(X, 'names', names, ...)         % static — design columns
%
% R = compute_rsm(fmri_data_obj, ...)                 % @fmri_data method
%
% -------------------------------------------------------------------------
% Author & copyright
% -------------------------------------------------------------------------
%   Copyright (C) 2026 Michael Sun.
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.

classdef rsm

    properties

        dat               = []
        is_dissimilarity  = false
        metric            = 'none'
        labels            = {}
        metadata_table    = table.empty
        groupings         = struct()
        level             = 'none'
        replicate_table   = table.empty
        whitened          = struct('level','none','method','none','shrinkage',[])
        source            = 'custom'
        history           = {}
        additional_info   = struct()

    end

    methods

        function obj = rsm(dat, varargin)
            % Construct an rsm from a k x k or k x k x N matrix plus optional
            % name-value metadata.
            %
            % Empty form:
            %   R = rsm();
            %
            % Bare-matrix form (defaults filled, labels auto-named):
            %   R = rsm(dat);
            %
            % Full form:
            %   R = rsm(dat, 'labels', {...}, 'metric', 'correlation', ...)

            if nargin == 0 || isempty(dat)
                return
            end

            if isnumeric(dat) || islogical(dat)
                obj.dat = double(dat);
            elseif isa(dat, 'rsm')
                % Copy-construct from another rsm; allow overrides via varargin.
                obj = dat;
            else
                error('rsm:badInput', ...
                    'First argument must be numeric, logical, or an rsm.');
            end

            % Parse name-value overrides
            p = inputParser;
            p.KeepUnmatched = false;
            p.CaseSensitive = false;
            p.addParameter('is_dissimilarity', obj.is_dissimilarity, @(x) islogical(x) || isnumeric(x));
            p.addParameter('metric',           obj.metric,           @(x) ischar(x) || isstring(x));
            p.addParameter('labels',           obj.labels,           @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
            p.addParameter('metadata_table',   obj.metadata_table,   @(x) istable(x) || isempty(x));
            p.addParameter('groupings',        obj.groupings,        @isstruct);
            p.addParameter('level',            obj.level,            @(x) ischar(x) || isstring(x));
            p.addParameter('replicate_table',  obj.replicate_table,  @(x) istable(x) || isempty(x));
            p.addParameter('whitened',         obj.whitened,         @isstruct);
            p.addParameter('source',           obj.source,           @(x) ischar(x) || isstring(x));
            p.addParameter('history',          obj.history,          @(x) iscell(x) || isempty(x));
            p.addParameter('additional_info',  obj.additional_info,  @isstruct);
            p.parse(varargin{:});

            obj.is_dissimilarity = logical(p.Results.is_dissimilarity);
            obj.metric           = char(p.Results.metric);
            obj.labels           = local_to_cellstr(p.Results.labels);
            obj.metadata_table   = p.Results.metadata_table;
            obj.groupings        = p.Results.groupings;
            obj.level            = char(p.Results.level);
            obj.replicate_table  = p.Results.replicate_table;
            obj.whitened         = p.Results.whitened;
            obj.source           = char(p.Results.source);
            obj.history          = p.Results.history;
            obj.additional_info  = p.Results.additional_info;

            % Auto-fill labels if missing
            k = size(obj.dat, 1);
            if isempty(obj.labels) && k > 0
                obj.labels = arrayfun(@(i) sprintf('cond_%d', i), 1:k, 'UniformOutput', false)';
            end

            % Run validators
            validate_metadata(obj);

            % Append a history entry
            obj.history{end+1} = sprintf('%s: constructed (%dx%dx%d, metric=%s, level=%s)', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), size(obj.dat,1), size(obj.dat,2), size(obj.dat,3), ...
                obj.metric, obj.level);

        end

    end

    methods (Static)
        % Static constructors live in separate files in @rsm/:
        %   from_categorical.m
        %   from_metadata_distance.m
        %   from_design.m
        %   from_table.m
        obj = from_categorical(varargin)
        obj = from_metadata_distance(varargin)
        obj = from_design(varargin)
        obj = from_table(varargin)
    end

end


% =========================================================================
% Local helpers (file-private)
% =========================================================================

function c = local_to_cellstr(x)
% Coerce label-like inputs to a column cellstr.
    if isempty(x)
        c = {};
    elseif iscellstr(x) %#ok<ISCLSTR>
        c = x(:);
    elseif isstring(x)
        c = cellstr(x(:));
    elseif ischar(x)
        c = {x};
    else
        error('rsm:badLabels', 'labels must be cellstr, string, or char.');
    end
end
