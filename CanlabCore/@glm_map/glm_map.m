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

        % --- Design specification (1st-level / event mode) ---------------
        design              % fmri_glm_design_matrix object (wrapped). Holds onsets, durations, event names/types, parametric modulators, basis set, and built X. Empty for pure direct/group designs.

        % --- Design specification (shared / direct mode) ----------------
        level = 2;          % Analysis level: 1 = first-level (within-run), 2 = second-level (group)
        is_timeseries = false;  % Logical; true if data are a within-run BOLD timeseries (enables AR error models in fit)
        contrasts = [];     % [regressors x contrasts] contrast matrix C
        contrast_names = {};    % Cell array of contrast names, one per column of C

        % --- Fitted result maps (populated by fit) ----------------------
        betas               % statistic_image (type 'Beta'), [voxels x regressors]
        t                   % statistic_image (type 'T') for regressors, [voxels x regressors]
        contrast_estimates  % statistic_image (type 'Contrast'), [voxels x contrasts]
        contrast_t          % statistic_image (type 'T') for contrasts, [voxels x contrasts]
        sigma               % fmri_data object: residual standard deviation per voxel
        dfe                 % Error degrees of freedom for the fit
        residuals           % fmri_data object: residuals [observations x voxels] (optional; only if requested)

        % --- Design diagnostics (populated by fit/diagnostics) ----------
        vif                 % Row vector of variance inflation factors, one per regressor (from VIF.m)
        contrast_vif        % Vector of contrast variance inflation factors, one per contrast (from cVIF.m)
        leverages           % Per-observation leverage values, diag(X*pinv(X))
        condition_number    % Condition number of the design matrix X
        rank_deficient = false;  % Logical; true if rank(X) < number of regressors
        collinearity_report     % Struct with redundant/duplicate row checks, near-collinear pairs, and centering flags
        warnings = {};      % Cell array of warning messages accumulated during build/fit

        % --- Provenance / metadata --------------------------------------
        analysis_name = '';     % Short descriptive name for this analysis
        fit_parameters = struct();  % Struct recording options used in fit (robust, AR order, threshold, grandmeanscale, ...)
        notes = '';         % Free-text notes
        history = {};       % Cell array of one-line provenance strings

        % --- Stored backing for direct-mode design ----------------------
        % (Use the Dependent .X / .regressor_names accessors instead of
        %  touching these directly.)
        Xdirect = [];           % [obs x regressors] design matrix for direct/group mode (no event model)
        regressor_names_direct = {};  % Cell array of regressor names for direct mode

        property_descriptions = { ...
            'design: fmri_glm_design_matrix object holding the 1st-level event design (onsets, durations, names, basis set, built X)' ...
            'level: 1 = first-level (within-run), 2 = second-level (group)' ...
            'is_timeseries: logical, true if data are a within-run BOLD timeseries (enables AR error models)' ...
            'contrasts: [regressors x contrasts] contrast matrix C' ...
            'contrast_names: cell array of contrast names' ...
            'betas/t: statistic_image maps of regression coefficients and their t-statistics, [voxels x regressors]' ...
            'contrast_estimates/contrast_t: statistic_image maps for linear contrasts, [voxels x contrasts]' ...
            'sigma: fmri_data object with residual standard deviation per voxel' ...
            'dfe: error degrees of freedom' ...
            'residuals: fmri_data object with residuals (optional)' ...
            'vif/contrast_vif: variance inflation factors for regressors and contrasts' ...
            'leverages: per-observation leverage values' ...
            'condition_number/rank_deficient: design matrix conditioning diagnostics' ...
            'collinearity_report: struct of redundant-row and near-collinearity checks' ...
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

            % If first argument is an fmri_glm_design_matrix, wrap it as a
            % 1st-level (event) design and consume that argument.
            if ~isempty(varargin) && isa(varargin{1}, 'fmri_glm_design_matrix')
                obj.design = varargin{1};
                obj.level = 1;
                varargin(1) = [];
            end

            % Names of stored (settable) properties for generic assignment
            stored_names = properties('glm_map');   % stored props only (Dependent excluded by properties())

            % Names of Dependent properties that have setters
            settable_dependent = {'TR', 'X', 'regressor_names'};

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
        % disp: concise object summary
        % =================================================================
        function disp(obj)

            fprintf('  glm_map object\n');
            fprintf('  %s\n', repmat('-', 1, 60));

            if ~isempty(obj.analysis_name)
                fprintf('  analysis_name : %s\n', obj.analysis_name);
            end

            switch obj.level
                case 1, levelstr = '1 (first-level / within-run)';
                case 2, levelstr = '2 (second-level / group)';
                otherwise, levelstr = num2str(obj.level);
            end
            fprintf('  level         : %s\n', levelstr);
            fprintf('  is_timeseries : %d\n', obj.is_timeseries);

            if ~isempty(obj.design)
                fprintf('  design        : fmri_glm_design_matrix (event mode), TR = %g\n', obj.design.TR);
            else
                fprintf('  design        : (direct/group mode, no event design)\n');
            end

            fprintf('  X             : %d images x %d regressors\n', obj.num_images, obj.num_regressors);
            fprintf('  contrasts     : %d\n', obj.num_contrasts);

            if obj.is_fitted
                fprintf('  fitted        : YES  (betas, t%s populated; dfe = %s)\n', ...
                    local_tf(~isempty(obj.contrast_estimates), ', contrasts'), num2str(obj.dfe));
            else
                fprintf('  fitted        : no   (run fit(obj, fmri_data_obj))\n');
            end

            if ~isempty(obj.vif)
                fprintf('  max VIF       : %.2f\n', max(obj.vif));
            end

            if ~isempty(obj.warnings)
                fprintf('  warnings      : %d (see obj.warnings)\n', numel(obj.warnings));
            end

            fprintf('  %s\n', repmat('-', 1, 60));
            fprintf('  methods(glm_map) for a list of operations.\n\n');

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


function s = local_tf(tf, str)
% Return str if tf is true, else ''.
if tf, s = str; else, s = ''; end
end % local_tf
