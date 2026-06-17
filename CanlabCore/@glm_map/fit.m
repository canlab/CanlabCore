function obj = fit(obj, data, varargin)
% Fit a mass-univariate GLM and populate result maps in a glm_map object.
%
% Scikit-learn-style fit: builds the design (if event/1st-level), runs
% fmri_data.regress on the data, and unpacks the outputs into the object's
% statistic_image result maps (betas, t, contrast_estimates, contrast_t),
% along with sigma, dfe, and the design diagnostics.
%
% :Usage:
% ::
%
%     obj = fit(obj, data, varargin)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object with a design specified (either a .design event
%        model or a direct .X matrix).
%
%   **data:**
%        An fmri_data object whose images correspond to the rows of X
%        (timeseries for level 1, contrast/subject images for level 2).
%
% :Optional Inputs:
%
%   **'robust':**
%        Use robust (bisquare) regression. Passed through to fmri_data.regress.
%
%   **'AR', [order]:**
%        Use an autoregressive error model of the given order. Only valid
%        when obj.is_timeseries is true.
%
%   **'residuals':**
%        Also store residuals in obj.residuals.
%
%   **'doverbose' / 'noverbose':**
%        Toggle verbose output (default true).
%
% :Outputs:
%
%   **obj:**
%        The glm_map object with result maps and diagnostics populated, and
%        is_fitted == true.
%
% :Examples:
% ::
%
%     dat = load_image_set('emotionreg');
%     g   = glm_map('X', [ones(30,1) zscore((1:30)')], 'level', 2);
%     g   = fit(g, dat);
%
% :See also:
%   - fmri_data.regress, build_design, diagnostics
%
% ..
%    Programmers' notes:
%    SCAFFOLD - not yet implemented. Planned behavior:
%      1. If event/1st-level and X not built, call build_design(obj).
%      2. Attach obj.X to data.X (and obj.contrasts/contrast_names).
%      3. If obj.is_timeseries and 'AR' requested, pass 'AR' to regress;
%         otherwise error if 'AR' requested on non-timeseries data.
%      4. out = regress(data, ...): unpack out.b -> obj.betas, out.t -> obj.t,
%         out.contrast_images -> obj.contrast_estimates, out.con_t ->
%         obj.contrast_t, out.sigma -> obj.sigma, out.df -> obj.dfe.
%      5. Populate obj.vif/contrast_vif/leverages/warnings from out.diagnostics
%         (or call diagnostics(obj)).
%      6. Record options in obj.fit_parameters and append to obj.history.
% ..

error('glm_map:NotImplemented', ...
    'fit() is scaffolded but not yet implemented. Planned for Phase 2 (wraps fmri_data.regress).');

end % fit
