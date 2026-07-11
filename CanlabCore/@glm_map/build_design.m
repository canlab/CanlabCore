function obj = build_design(obj, varargin)
% Build the design matrix X for an event/1st-level glm_map from onsets.
%
% Delegates to the wrapped fmri_glm_design_matrix object's build method,
% which convolves onsets/durations with the basis set and assembles the
% design matrix (interest, covariate, and baseline partitions) into
% design.xX.X with column names in design.xX.name. After this call, the
% Dependent obj.X and obj.regressor_names read through to the built design.
%
% :Usage:
% ::
%
%     obj = build_design(obj, varargin)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object with a non-empty .design (fmri_glm_design_matrix)
%        that has conditions/onsets assigned.
%
% :Optional Inputs:
%
%   **'doplot' / 'plot':**
%        Plot the resulting design via fmri_glm_design_matrix.plot (default off).
%
% :Outputs:
%
%   **obj:**
%        The glm_map object with the wrapped design built (design.xX.X populated).
%
% :Examples:
% ::
%
%     TR = 2; nscan = 200;
%     ons = {[10 40 70 100]', [25 55 85 115]'};   % 2 conditions, 1 session
%     d = fmri_glm_design_matrix(TR, 'nscan', nscan, 'units', 'secs', ...
%             'onsets', ons, 'condition_names', {'A' 'B'});
%     g = glm_map(d);
%     g = build_design(g);
%     size(g.X)            % built design matrix
%     g.regressor_names
%
% :See also:
%   - fmri_glm_design_matrix, fmri_glm_design_matrix.build, onsets2fmridesign
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation (delegates to fmri_glm_design_matrix.build).
% ..

% -------------------------------------------------------------------------
% Parse options
% -------------------------------------------------------------------------
doplot = any(strcmpi(varargin, 'doplot')) || any(strcmpi(varargin, 'plot'));

% -------------------------------------------------------------------------
% Validate
% -------------------------------------------------------------------------
if isempty(obj.design)
    error('glm_map:NoDesign', ...
        'build_design requires an fmri_glm_design_matrix in obj.design (event/1st-level mode).');
end

if isempty(obj.design.Sess) || isempty(obj.design.Sess(1).U) || isempty(obj.design.Sess(1).U(1).name)
    error('glm_map:NoOnsets', ...
        ['The wrapped design has no conditions/onsets assigned. Add them first, e.g. ' ...
         'fmri_glm_design_matrix(TR, ''nscan'', nscan, ''onsets'', ons, ''condition_names'', names).']);
end

% -------------------------------------------------------------------------
% Build
% -------------------------------------------------------------------------
obj.design = build(obj.design);

obj.level = 1;   % building implies a first-level/event model

if doplot
    plot(obj.design);
end

obj.history{end + 1} = sprintf('build_design: built design matrix [%d x %d] from onsets', ...
    size(obj.X, 1), size(obj.X, 2));

end % build_design
