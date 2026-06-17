function obj = build_design(obj, varargin)
% Build the design matrix X for an event/1st-level glm_map from onsets.
%
% Delegates to the wrapped fmri_glm_design_matrix object (and ultimately
% onsets2fmridesign) to convolve onsets/durations with the basis set and
% assemble X, regressor names, and any nuisance covariates. After this call,
% the Dependent obj.X and obj.regressor_names read through to the built design.
%
% :Usage:
% ::
%
%     obj = build_design(obj, varargin)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object with a non-empty .design (fmri_glm_design_matrix).
%
% :Optional Inputs:
%
%   **'doplot' / 'noplot':**
%        Plot the resulting design matrix (default false).
%
% :Outputs:
%
%   **obj:**
%        The glm_map object with the wrapped design built (design.xX.X populated).
%
% :See also:
%   - fmri_glm_design_matrix, fmri_glm_design_matrix.build, onsets2fmridesign
%
% ..
%    Programmers' notes:
%    SCAFFOLD - not yet implemented. Planned behavior:
%      1. Require obj.level == 1 and ~isempty(obj.design).
%      2. obj.design = build(obj.design);   % existing method
%      3. Sync obj.contrasts/contrast_names defaults if requested.
%      4. Append to obj.history.
% ..

if isempty(obj.design)
    error('glm_map:NoDesign', ...
        'build_design requires an fmri_glm_design_matrix in obj.design (event/1st-level mode).');
end

error('glm_map:NotImplemented', ...
    'build_design() is scaffolded but not yet implemented. Planned for Phase 3 (wraps fmri_glm_design_matrix.build).');

end % build_design
