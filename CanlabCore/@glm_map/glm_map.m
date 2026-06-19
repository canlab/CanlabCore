% glm_map Object class for mass-univariate GLM / multiple-regression fits.
%
% A scikit-learn-style estimator object that bundles, in a single container:
%   (1) the *design specification* (onsets, durations, event names/types,
%       parametric modulators, basis set, covariates, design matrix X, and
%       contrasts) -- wrapping an fmri_glm_design_matrix object,
%   (2) the *fitted result maps* (betas, t, contrast estimates, contrast t),
%       stored as statistic_image objects, plus sigma and residuals, and
%   (3) the *design diagnostics* (variance inflation factors, contrast VIFs,
%       leverage, condition number, rank, collinearity/redundancy checks).
%
% glm_map holds the outputs of a mass-univariate GLM fit produced by
% fmri_data.regress. The workflow mirrors scikit-learn estimators:
%
%   g = glm_map(...)            % construct: set design + fit options (no data)
%   g = fit(g, fmri_data_obj)   % build design (if needed) and run regression
%   table(g); montage(g);       % inspect results
%
% - This is a standalone container class (composition, not an image_vector
%   subclass): the voxel data live inside the contained statistic_image /
%   fmri_data objects (betas, t, contrast_estimates, ...).
% - It supports two design modes:
%     * level 1 (event / 1st-level): an fmri_glm_design_matrix in .design
%       builds X by HRF convolution of onsets. Mark .is_timeseries = true
%       for within-run BOLD data so autoregressive ('AR') models are valid.
%     * level 2 (direct / group): a design matrix X is supplied directly.
% - Several convenience attributes (X, regressor_names, onsets, durations,
%   condition_names, TR, num_*) are implemented as true MATLAB *Dependent*
%   properties that read through to the wrapped design object, so there is a
%   single source of truth and no risk of stale duplicates.
%
% :Usage:
% ::
%
%     g = glm_map(varargin)
%     g = glm_map(fmri_glm_design_matrix_obj)         % wrap a 1st-level design
%     g = glm_map('X', X, 'level', 2)                 % direct/group design
%     g = glm_map('fieldname', value, ...)            % set any stored property
%     g = glm_map(regress_out_struct)                 % re-cast a regress() output struct
%
% The object property names mirror the fields of the results structure
% returned by fmri_data.regress (the variable "out" in that method). Related
% options are grouped into nested structs (.input_parameters,
% .input_image_metadata, .diagnostics), and per-voxel df / sigma are kept as
% fmri_data objects in .df and .sigma. For backward compatibility the
% historical out-struct field names are available as aliases that read/write
% the canonical properties:
%
%     .b -> .betas    .contrast_images -> .contrast_estimates    .con_t -> .contrast_t
%     .resid -> .residuals    .variable_names -> .regressor_names    .C -> .contrasts
%
% :Inputs:
%
%   **(optional first arg) fmri_glm_design_matrix object:**
%        If the first argument is an fmri_glm_design_matrix, it is stored in
%        .design and the object is marked as a 1st-level (event) model.
%
%   **'fieldname', value pairs:**
%        Any settable property of glm_map followed by a value. See
%        properties(glm_map) and the property_descriptions field.
%
% :Optional Inputs:
%
%   **'X', [obs x regressors] matrix:**
%        A pre-built design matrix for direct/group (level-2) analysis.
%
%   **'level', 1 | 2:**
%        Analysis level. 1 = first-level (within-run); 2 = second-level
%        (group). Default 2 for direct designs, 1 when a design object is given.
%
%   **'is_timeseries', [logical]:**
%        True if the data to be fit are a within-run BOLD timeseries; enables
%        autoregressive (AR) error models in fit. Default false.
%
%   **'contrasts', [regressors x contrasts] matrix:**
%        Contrast matrix C (rows must match number of regressors in X).
%
%   **'contrast_names', {cellstr}:**
%        Names for each contrast (columns of C).
%
% :Outputs:
%
%   **g:**
%        A glm_map object.
%
% :Examples:
% ::
%
%     % ---- Direct / group (2nd-level) design ----
%     dat = load_image_set('emotionreg');          % 30 contrast images
%     X = [ones(30,1) zscore((1:30)')];            % intercept + a covariate
%     g = glm_map('X', X, 'level', 2, ...
%                 'regressor_names', {'intercept' 'cov'});
%     g = fit(g, dat);                             % runs fmri_data.regress
%     diagnostics(g);                              % VIFs, leverage, etc.
%     table(g); montage(g, 't');
%
%     % ---- Event / 1st-level design wrapping fmri_glm_design_matrix ----
%     d = fmri_glm_design_matrix(TR, 'nscan', nscan, 'units', 'secs', ...
%             'onsets', ons, 'condition_names', names);
%     g = glm_map(d);                              % level 1, event mode
%     g.is_timeseries = true;                      % enable AR models in fit
%     g = build_design(g);                         % onsets -> X via convolution
%     g.onsets                                     % read-through to .design
%
%     % ---- Import from an SPM (SPM12/SPM25) first-level model ----
%     g = import_SPM(glm_map, '/path/to/SPM.mat');
%
% :See also:
%   - fmri_data, fmri_data.regress
%   - statistic_image
%   - fmri_glm_design_matrix
%   - VIF, cVIF
%
% ..
%    Author and copyright information:
%
%    Copyright (C) 2026 Tor Wager
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
%    Programmers' notes:
%    2026 - Initial scaffold. Constructor, Dependent accessors, and disp are
%    functional; fit/build_design/diagnostics/import_SPM/display/table are
%    stubbed for later implementation phases.
% ..

classdef glm_map

    % ---------------------------------------------------------------------
    % Stored properties
    % ---------------------------------------------------------------------
    properties

        % --- Provenance / metadata --------------------------------------
        analysis_name = '';     % Short descriptive name for this analysis (mirrors regress out.analysis_name)

        % --- Design specification (1st-level / event mode) ---------------
        design              % fmri_glm_design_matrix object (wrapped). Holds onsets, durations, event names/types, parametric modulators, basis set, and built X. Empty for pure direct/group designs.

        % --- Design specification (shared / direct mode) ----------------
        level = 2;          % Analysis level: 1 = first-level (within-run), 2 = second-level (group)
        is_timeseries = false;  % Logical; true if data are a within-run BOLD timeseries (enables AR error models in fit)
        contrasts = [];     % [regressors x contrasts] contrast matrix C (read/written as out.C via the .C alias)
        contrast_names = {};    % Cell array of contrast names, one per column of C (mirrors out.contrast_names)
        contrast_summary_table = table();  % Table summarizing contrasts over conditions (mirrors out.contrast_summary_table)

        % --- Nested input/option structs (mirror fmri_data.regress out) --
        input_parameters = struct();      % Struct: options used by the fit (brain_is_predictor, do_robust, grandmeanscale, do_intercept, do_resid, initial_statistical_threshold, ...). Mirrors out.input_parameters.
        input_image_metadata = struct();  % Struct: provenance of the fitted images (source_notes, history, image_names, fullpath). Mirrors out.input_image_metadata.

        % --- Fitted result maps (populated by fit / ingested from regress) -
        betas               % statistic_image (type 'Beta'), [voxels x regressors].   Alias: .b
        t                   % statistic_image (type 'T') for regressors, [voxels x regressors]
        contrast_estimates  % statistic_image (type 'Contrast'), [voxels x contrasts]. Alias: .contrast_images
        contrast_t          % statistic_image (type 'T') for contrasts, [voxels x contrasts]. Alias: .con_t
        df                  % fmri_data object: per-voxel error degrees of freedom (mirrors out.df)
        sigma               % fmri_data object: residual standard deviation per voxel (mirrors out.sigma)
        residuals           % fmri_data object: residuals (optional; only if requested). Alias: .resid
        dfe                 % Scalar error degrees of freedom for the fit (median of df, convenience summary)

        % --- Design diagnostics (nested struct; populated by fit/diagnostics) -
        % Fields: Variance_inflation_factors, Leverages (same names as
        % fmri_data.regress out.diagnostics), plus
        % Contrast_variance_inflation_factors, Cooks_distance (per-observation
        % influence; populated when residuals are available), condition_number,
        % rank_deficient, collinearity_report, vif_threshold.
        diagnostics = struct();

        warnings = {};      % Cell array of warning messages accumulated during build/fit (mirrors out.warnings)

        % --- Provenance / metadata --------------------------------------
        fit_parameters = struct();  % Struct recording options used in fit (robust, AR order, threshold, ...)
        notes = '';         % Free-text notes
        history = {};       % Cell array of one-line provenance strings

        % --- Stored backing for direct-mode design ----------------------
        % (Use the Dependent .X / .regressor_names accessors instead of
        %  touching these directly.)
        Xdirect = [];           % [obs x regressors] design matrix for direct/group mode (no event model)
        regressor_names_direct = {};  % Cell array of regressor names for direct mode

        property_descriptions = { ...
            'analysis_name: short descriptive name (out.analysis_name)' ...
            'design: fmri_glm_design_matrix object holding the 1st-level event design (onsets, durations, names, basis set, built X)' ...
            'level: 1 = first-level (within-run), 2 = second-level (group)' ...
            'is_timeseries: logical, true if data are a within-run BOLD timeseries (enables AR error models)' ...
            'contrasts (.C): [regressors x contrasts] contrast matrix' ...
            'contrast_names: cell array of contrast names' ...
            'contrast_summary_table: table of contrast weights over conditions (out.contrast_summary_table)' ...
            'input_parameters: struct of fit options (out.input_parameters)' ...
            'input_image_metadata: struct of input-image provenance (out.input_image_metadata)' ...
            'betas (.b) / t: statistic_image maps of regression coefficients and their t-statistics, [voxels x regressors]' ...
            'contrast_estimates (.contrast_images) / contrast_t (.con_t): statistic_image maps for linear contrasts, [voxels x contrasts]' ...
            'df: fmri_data object with per-voxel error degrees of freedom (out.df)' ...
            'sigma: fmri_data object with residual standard deviation per voxel (out.sigma)' ...
            'residuals (.resid): fmri_data object with residuals (optional)' ...
            'dfe: scalar error degrees of freedom (median of df)' ...
            'diagnostics: struct with vif, contrast_vif, leverages, condition_number, rank_deficient, collinearity_report' ...
            'warnings: cell array of warnings from build/fit' ...
            'analysis_name/fit_parameters/notes/history: provenance and metadata' ...
            };

    end % stored properties


    % ---------------------------------------------------------------------
    % Dependent properties (true MATLAB Dependent; computed accessors)
    % ---------------------------------------------------------------------
    properties (Dependent)

        TR                  % Repetition time (s). Reads/writes design.TR for event models.
        X                   % [obs x regressors] design matrix. Reads built X from .design (event) or .Xdirect (direct).
        regressor_names     % Cell array of regressor (design column) names
        onsets              % Cell array of event onsets (read-through to .design; 1st-level only)
        durations           % Cell array of event durations (read-through to .design; 1st-level only)
        condition_names     % Cell array of condition/event names (read-through to .design; 1st-level only)
        num_images          % Number of images/observations (rows of X)
        num_regressors      % Number of regressors (columns of X)
        num_contrasts       % Number of contrasts (columns of C)
        is_fitted           % Logical; true once result maps (.betas) are populated

        % --- Aliases to the fmri_data.regress out-struct field names -----
        % These read/write the canonical properties above, so that an object
        % returned by fmri_data.regress supports the historical struct-style
        % field access (out.b, out.con_t, ...) unchanged.
        b                   % Alias for betas               (out.b)
        contrast_images     % Alias for contrast_estimates  (out.contrast_images)
        con_t               % Alias for contrast_t          (out.con_t)
        resid               % Alias for residuals           (out.resid)
        variable_names      % Alias for regressor_names     (out.variable_names)
        C                   % Alias for contrasts           (out.C)

    end % dependent properties


    methods

        % =================================================================
        % Constructor
        % =================================================================
        function obj = glm_map(varargin)

            % Empty object: return defaults
            if nargin == 0
                return
            end

            % If the first argument is a struct, treat it as a regression
            % results structure (the historical fmri_data.regress output) and
            % re-cast it as a glm_map object. Any remaining 'field', value
            % pairs are then applied as overrides.
            if ~isempty(varargin) && isstruct(varargin{1})
                obj = local_from_regress_struct(obj, varargin{1});
                varargin(1) = [];

            % If first argument is an fmri_glm_design_matrix, wrap it as a
            % 1st-level (event) design and consume that argument.
            elseif ~isempty(varargin) && isa(varargin{1}, 'fmri_glm_design_matrix')
                obj.design = varargin{1};
                obj.level = 1;
                varargin(1) = [];
            end

            % Names of stored (settable) properties for generic assignment
            stored_names = properties('glm_map');   % includes public Dependent props

            % Names of Dependent properties that have setters
            settable_dependent = {'TR', 'X', 'regressor_names', ...
                'b', 'contrast_images', 'con_t', 'resid', 'variable_names', 'C'};

            for i = 1:length(varargin)

                if ~ischar(varargin{i}), continue, end

                fieldname = varargin{i};

                % Map a couple of friendly aliases to backing storage
                switch fieldname
                    case {'names', 'variable_names'}
                        fieldname = 'regressor_names';
                end

                if any(strcmp(fieldname, stored_names)) || any(strcmp(fieldname, settable_dependent))

                    obj.(fieldname) = varargin{i + 1};

                    % Blank the consumed value so a trailing char value is
                    % not re-interpreted as a keyword on the next iteration
                    if ischar(varargin{i + 1})
                        varargin{i + 1} = [];
                    end

                else
                    warning('glm_map:UnknownField', 'Unknown glm_map field: %s', fieldname);
                end

            end % parse inputs

        end % constructor


        % =================================================================
        % Dependent property GET accessors
        % =================================================================
        function val = get.TR(obj)
            if ~isempty(obj.design)
                val = obj.design.TR;
            else
                val = [];
            end
        end

        function val = get.X(obj)
            % Prefer a built design matrix from the wrapped design object;
            % fall back to the direct-mode matrix.
            val = [];
            if ~isempty(obj.design) && isstruct(obj.design.xX) && isscalar(obj.design.xX) ...
                    && isfield(obj.design.xX, 'X') && ~isempty(obj.design.xX.X)
                val = obj.design.xX.X;
            elseif ~isempty(obj.Xdirect)
                val = obj.Xdirect;
            end
        end

        function val = get.regressor_names(obj)
            val = {};
            if ~isempty(obj.design) && isstruct(obj.design.xX) && isscalar(obj.design.xX) ...
                    && isfield(obj.design.xX, 'name') && ~isempty(obj.design.xX.name)
                val = obj.design.xX.name;
            elseif ~isempty(obj.regressor_names_direct)
                val = obj.regressor_names_direct;
            end
        end

        function val = get.onsets(obj)
            val = local_collect_U_field(obj.design, 'ons');
        end

        function val = get.durations(obj)
            val = local_collect_U_field(obj.design, 'dur');
        end

        function val = get.condition_names(obj)
            val = local_collect_U_field(obj.design, 'name');
        end

        function val = get.num_images(obj)
            val = size(obj.X, 1);
        end

        function val = get.num_regressors(obj)
            val = size(obj.X, 2);
        end

        function val = get.num_contrasts(obj)
            val = size(obj.contrasts, 2);
        end

        function val = get.is_fitted(obj)
            val = ~isempty(obj.betas);
        end


        % =================================================================
        % Dependent property SET accessors
        % =================================================================
        function obj = set.TR(obj, val)
            if isempty(obj.design)
                error('glm_map:NoDesign', ...
                    'Cannot set TR: no fmri_glm_design_matrix in .design. Create one first (event/1st-level mode).');
            end
            obj.design.TR = val;
        end

        function obj = set.X(obj, val)
            % Setting X targets the direct/group design backing store.
            if ~isempty(obj.design) && isstruct(obj.design.xX) && isscalar(obj.design.xX) ...
                    && isfield(obj.design.xX, 'X') && ~isempty(obj.design.xX.X)
                warning('glm_map:DesignPresent', ...
                    ['This glm_map has a built event design in .design; setting .X stores a direct ' ...
                     'matrix in .Xdirect that the event design will shadow. Use build_design instead.']);
            end
            obj.Xdirect = val;
        end

        function obj = set.regressor_names(obj, val)
            if ~iscell(val), val = cellstr(val); end
            obj.regressor_names_direct = val;
        end


        % =================================================================
        % Aliases to the fmri_data.regress out-struct field names
        % (read/write the canonical properties)
        % =================================================================
        function val = get.b(obj),               val = obj.betas;              end
        function obj = set.b(obj, val),           obj.betas = val;              end

        function val = get.contrast_images(obj),  val = obj.contrast_estimates; end
        function obj = set.contrast_images(obj, val), obj.contrast_estimates = val; end

        function val = get.con_t(obj),            val = obj.contrast_t;         end
        function obj = set.con_t(obj, val),       obj.contrast_t = val;         end

        function val = get.resid(obj),            val = obj.residuals;          end
        function obj = set.resid(obj, val),       obj.residuals = val;          end

        function val = get.variable_names(obj),   val = obj.regressor_names;    end
        function obj = set.variable_names(obj, val), obj.regressor_names = val; end

        function val = get.C(obj),                val = obj.contrasts;          end
        function obj = set.C(obj, val),           obj.contrasts = val;          end


        % =================================================================
        % disp: object summary with full property listing
        % =================================================================
        function disp(obj)

            line = repmat('-', 1, 64);
            fprintf('  glm_map object\n  %s\n', line);

            % -------- One-line status header --------
            switch obj.level
                case 1, levelstr = '1 (first-level / within-run)';
                case 2, levelstr = '2 (second-level / group)';
                otherwise, levelstr = num2str(obj.level);
            end
            fprintf('  level %s | X: %d images x %d regressors | %d contrast(s) | fitted: %s\n', ...
                levelstr, obj.num_images, obj.num_regressors, obj.num_contrasts, ...
                local_tf(obj.is_fitted, 'YES', 'no'));
            fprintf('  %s\n', line);

            % -------- Full property listing --------
            % Curated order; nested structs are expanded one level so their
            % fields are visible at a glance. (See properties(obj) for the
            % full set, including the out-struct aliases b/con_t/.../C.)
            props = { ...
                'analysis_name', 'design', 'level', 'is_timeseries', ...
                'contrasts', 'contrast_names', 'contrast_summary_table', ...
                'input_parameters', 'input_image_metadata', ...
                'betas', 't', 'contrast_estimates', 'contrast_t', ...
                'df', 'sigma', 'residuals', 'dfe', ...
                'diagnostics', 'warnings', ...
                'fit_parameters', 'notes', 'history'};

            for i = 1:numel(props)
                p = props{i};
                fprintf('  %-22s : %s\n', p, local_summarize(obj.(p)));

                % Expand the three nested option/diagnostic structs one level
                if any(strcmp(p, {'input_parameters', 'input_image_metadata', 'diagnostics'})) ...
                        && isstruct(obj.(p)) && ~isempty(fieldnames(obj.(p)))
                    fn = fieldnames(obj.(p));
                    for j = 1:numel(fn)
                        fprintf('  %-22s     .%-18s %s\n', '', fn{j}, ...
                            local_summarize(obj.(p).(fn{j})));
                    end
                end
            end

            fprintf('  %s\n', line);
            fprintf('  methods(glm_map) for operations; properties(glm_map) for all fields.\n\n');

        end % disp

    end % methods

end % classdef


% =====================================================================
% Local helper functions
% =====================================================================
function out = local_collect_U_field(design, fieldname)
% Read-through accessor: collect a field (ons/dur/name) from the nested
% design.Sess(s).U(i) structure into a flat cell array. Returns {} if the
% wrapped design has no session/onset structure yet.

out = {};
if isempty(design) || ~isprop(design, 'Sess') || isempty(design.Sess)
    return
end

for s = 1:numel(design.Sess)
    if ~isfield(design.Sess(s), 'U') || isempty(design.Sess(s).U)
        continue
    end
    for u = 1:numel(design.Sess(s).U)
        if isfield(design.Sess(s).U(u), fieldname)
            out{end + 1} = design.Sess(s).U(u).(fieldname); %#ok<AGROW>
        end
    end
end

end % local_collect_U_field


function s = local_tf(tf, a, b)
% Inline if: return a if tf is true, else b ('' if b omitted).
if nargin < 3, b = ''; end
if tf, s = a; else, s = b; end
end % local_tf


function s = local_summarize(v)
% Compact one-line description of any property value, for disp().

if isempty(v)
    s = '[]';
    return
end

if ischar(v)
    if size(v, 1) > 1
        s = sprintf('[%dx%d char]', size(v, 1), size(v, 2));   % char matrix
    else
        s = v;
        if numel(s) > 52, s = [s(1:49) '...']; end
    end
    return
end

if islogical(v) && isscalar(v)
    s = local_tf(v, 'true', 'false');
    return
end

if isnumeric(v) && isscalar(v)
    s = num2str(v);
    return
end

if isa(v, 'statistic_image') || isa(v, 'fmri_data')
    try
        s = sprintf('%s [%d voxels x %d images]', class(v), size(v.dat, 1), size(v.dat, 2));
    catch
        s = class(v);
    end
    return
end

if isa(v, 'fmri_glm_design_matrix')
    try
        s = sprintf('fmri_glm_design_matrix (TR = %g)', v.TR);
    catch
        s = 'fmri_glm_design_matrix';
    end
    return
end

if istable(v)
    s = sprintf('table [%d x %d]', size(v, 1), size(v, 2));
    return
end

if isstruct(v)
    fn = fieldnames(v);
    if isempty(fn)
        s = 'struct (empty)';
    else
        s = sprintf('struct with %d field(s)', numel(fn));
    end
    return
end

if iscell(v)
    % Show contents if it is a short cellstr
    if numel(v) <= 6 && all(cellfun(@(x) ischar(x) || isstring(x), v))
        s = ['{' strjoin(cellfun(@char, v(:)', 'UniformOutput', false), ', ') '}'];
        if numel(s) > 52, s = sprintf('{%dx%d cell}', size(v, 1), size(v, 2)); end
    else
        s = sprintf('{%dx%d cell}', size(v, 1), size(v, 2));
    end
    return
end

if isnumeric(v)
    s = sprintf('[%dx%d %s]', size(v, 1), size(v, 2), class(v));
    return
end

s = sprintf('[%s]', class(v));

end % local_summarize


function obj = local_from_regress_struct(obj, S)
% Re-cast a regression results structure (the historical fmri_data.regress
% output) into a glm_map object. Maps the out-struct field names onto the
% canonical glm_map properties; tolerates missing fields.

obj.level = 2;   % regress output carries no level; default to group

% src field in struct  ->  dst glm_map property (or settable alias)
map = { ...
    'analysis_name',          'analysis_name'; ...
    'input_parameters',       'input_parameters'; ...
    'input_image_metadata',   'input_image_metadata'; ...
    'X',                      'X'; ...                 % set.X -> Xdirect
    'variable_names',         'regressor_names'; ...   % set.regressor_names -> regressor_names_direct
    'C',                      'contrasts'; ...
    'contrast_names',         'contrast_names'; ...
    'contrast_summary_table', 'contrast_summary_table'; ...
    'warnings',               'warnings'; ...
    'b',                      'betas'; ...
    't',                      't'; ...
    'df',                     'df'; ...
    'sigma',                  'sigma'; ...
    'resid',                  'residuals'; ...
    'contrast_images',        'contrast_estimates'; ...
    'con_t',                  'contrast_t'};

for i = 1:size(map, 1)
    src = map{i, 1};
    if ~isfield(S, src), continue, end
    v = S.(src);
    % Skip genuinely-empty primitives, but always ingest result-map objects
    % (statistic_image / fmri_data may overload isempty based on volInfo).
    if (isnumeric(v) || ischar(v) || iscell(v) || isstruct(v)) && isempty(v)
        continue
    end
    obj.(map{i, 2}) = v;
end

% Diagnostics: ingest as-is (regress already names the fields
% Variance_inflation_factors / Leverages, which glm_map reuses)
if isfield(S, 'diagnostics') && ~isempty(S.diagnostics)
    obj.diagnostics = S.diagnostics;
end

% Scalar error-df summary from the per-voxel df image
if isfield(S, 'df') && ~isempty(S.df) && isprop(S.df, 'dat') && ~isempty(S.df.dat)
    obj.dfe = double(median(S.df.dat(:)));
end

obj.history{end + 1} = 'constructed from fmri_data.regress() output structure';

end % local_from_regress_struct


