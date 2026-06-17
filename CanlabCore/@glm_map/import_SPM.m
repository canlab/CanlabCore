function obj = import_SPM(obj, SPM, varargin)
% Import a first-level design (and optional betas) from an SPM model.
%
% Populates a glm_map object from an SPM structure or an SPM.mat path,
% mapping SPM's design fields into the wrapped fmri_glm_design_matrix object.
% Compatible with SPM12 and SPM25 first-level models.
%
% :Usage:
% ::
%
%     obj = import_SPM(glm_map, SPM)            % SPM struct already loaded
%     obj = import_SPM(glm_map, '/path/SPM.mat')
%
% :Inputs:
%
%   **obj:**
%        A glm_map object (typically empty: glm_map).
%
%   **SPM:**
%        Either an SPM structure (as loaded from SPM.mat) or a char/string
%        path to an SPM.mat file.
%
% :Optional Inputs:
%
%   **'load_betas':**
%        Also load the estimated beta_*.nii images into obj.betas.
%
% :Outputs:
%
%   **obj:**
%        glm_map with .design populated from SPM, level set to 1, and
%        is_timeseries set true.
%
% :See also:
%   - fmri_glm_design_matrix, build_design
%
% ..
%    Programmers' notes:
%    SCAFFOLD - not yet implemented. Planned field mapping (SPM -> object):
%      SPM.Sess         -> obj.design.Sess   (onsets .ons, durations .dur,
%                          names .name, parametric mods .P, covariates .C)
%      SPM.xBF          -> obj.design.xBF    (basis set, T, T0, UNITS, Volterra)
%      SPM.xX           -> obj.design.xX     (design matrix .X, partitions, names)
%      SPM.xY.RT        -> obj.design.TR
%      SPM.nscan        -> obj.design.nscan
%      obj.level = 1; obj.is_timeseries = true;
%    SPM25 deltas vs SPM12 will be reconciled against a real SPM.mat at
%    implementation time.
% ..

error('glm_map:NotImplemented', ...
    'import_SPM() is scaffolded but not yet implemented. Planned for Phase 3.');

end % import_SPM
