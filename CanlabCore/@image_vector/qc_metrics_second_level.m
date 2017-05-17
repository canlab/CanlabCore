function [group_metrics individual_metrics values gwcsf gwcsfmean gwcsf_l2norm] = qc_metrics_second_level(obj, varargin)
%
% Quality metrics for a 2nd-level analysis (set of images from different subjects)
% The goal is to obtain measures that are simple and scale invariant, and
% so can be compared across image datasets
%
% :Usage:
% ::
%
%     [group_metrics individual_metrics] = qc_metrics_second_level(list inputs here, [optional inputs])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2016 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **obj:**
%        An fmri_data object with a set of images
%        Intended to be individual-participant beta or contrast images from
%        first-level analyses.
%
% :Optional Inputs:
%   **noverbose:**
%        Turn off text display reporting of warnings and metrics
%
%   **param2:**
%        description of param2
%
% :Outputs:
%
%   **group_metrics:**
%        Quality metrics indicating potential problems, as described below
%
%   **individual_metrics:**
%        Individual image-level values for quality metrics
%        These can be used to identify individuals with problematic images
%        or as covariates in analyses.
%
%   **values:**
%        Grand Mean Gray, white, and CSF values - one number per metric per
%        image
%
%   **gwcsf:**
%        Gray, white, and CSF image objects with values for all images
%
%   **gwcsfmean:**
%        Gray, white, and CSF image objects mean values across images
%
%   **gwcsfl2norm:**
%        L2 norm of Gray, white, and CSF components for each image
%
% :References:
%   None.
%
%
% --------------------------------------------------------------------------------
% Bright-edge bias:  "brightedge"
%       NOT IMPLEMENTED YET
%       High inter-image variability around edge of brain
%       median absolute deviation (MAD) in brain edge-mask / MAD in center
%       of brain.
%       Range: [0 Inf],  1 is good (homogeneous). High values may indicate
%       head movement artifacts
%
% Signal in ventricles: "csf_to_gm_signal_ratio"
%       In most contrast maps, there should be less non-zero signal in the CSF space (e.g., ventricles)
%       than in gray matter.  This measures the ratio of mean absolute signal in CSF / GM.
%       Range: [0 Inf]. Values < 1 indicate low ventricle signal, which is good.
%       Higher values are worse.
%
% Effect size in ventricles: "global_d_ventricles"
%       d = Mean / standard error (STD) of whole-ventricle average signal
%       Indicates a shift towards global activation or deactivation of the
%       ventricles in the group.
%       Range: [-Inf Inf], farther from zero is bad.
%
% Significant global activation of ventricles: "global_logp_ventricles"
%       -log(p-value) in ventricle global mean across subjects, divided by -log(0.05)
%       Any value > 1 indicates a significant global activation in the
%       ventricles.
%       Range: [0 Inf], higher is worse. Values > 1 are bad.
%
% Effect size in white matter: "global_d_wm"
%       d = Mean / standard error (STD) of whole-WM average signal
%       Indicates a shift towards global activation or deactivation of the
%       WM in the group. This is not completely redundant with CSF, because
%       some artifacts can affect one tissue compartment but not another.
%       WM values are, however, more likely to be influenced by true
%       signal.
%       Range: [-Inf Inf], farther from zero is bad.
%
% Significant global activation of WM: "global_logp_wm"
%       -log(p-value) in WM global mean across subjects, divided by -log(0.05)
%       Any value > 1 indicates a significant global activation in WM.
%       Range: [0 Inf], higher is worse. Values > 1 are bad.
%
% Non-spatially specific signal contamination: "r2_explained_by_csf"
%       Variance in individual differences in whole-gray-matter activity
%       Range: [-1 1], usually [0 1].  Farther from zero is bad.
%
% Scale inhomogeneity: "csf_scale_inhom"
%       The coefficient of variation of the L1-norm across CSF voxels for each participant.
%       The L1-norm is a measure of scale (across voxels, for each observation/image)
%       that is more robust than the variance.
%       High variability in these individual scale parameters is a sign of inhomogeneity.
%       The coefficient of variation is the standard dev. of this whole-image scale value divided by its mean.
%       This is a scale-invariant measure of inhomogeneity across images.
%       It also avoids the instability that occurs by dividing by values that may
%       be near zero, such as the mean.
%       Range: [0 Inf].  Lower is better.
%
% Scale inhomogeneity in gray matter: "gm_scale_inhom"
%       The coefficient of variation of the L1-norm across GM voxels for each participant.
%       The L1-norm is a measure of scale (across voxels, for each observation/image)
%       that is more robust than the variance.
%       High variability in these individual scale parameters is a sign of inhomogeneity.
%       The coefficient of variation is the standard dev. of this whole-image scale value divided by its mean.
%       This is a scale-invariant measure of inhomogeneity across images.
%       It also avoids the instability that occurs by dividing by values that may
%       be near zero, such as the mean.
%       Range: [0 Inf].  Lower is better, but GM values may also be influenced by large differences in true signal.
%
% :Examples:
% ::
%    % Load sample images
%    obj = load_image_set('emotionreg');
%
%    % Run QC metrics
%    [group_metrics individual_metrics] = qc_metrics_second_level(obj);
%
%    % Plot histograms of gray, white, CSF values
%    hist_han = histogram(obj, 'byimage', 'by_tissue_type');
%
% :See also:
%   - extract_gray_white_csf, histogram
%


% ..
%    DEFAULTS AND INPUTS
% ..


doverbose = 1;

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'noverbose', doverbose = 0;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


[values, components, gwcsf, gwcsf_l2norm] = extract_gray_white_csf(obj);

% Get mean signal in each tissue compartment for visualization, if desired
% -------------------------------------------------------------------------
gwcsfmean = mean(gwcsf{1}); 
gwcsfmean = replace_empty(gwcsfmean); 
for i = 2:3
    m = mean(gwcsf{i}); 
    m = replace_empty(m);
    gwcsfmean.dat(:, i) = m.dat;
end


% -------------------------------------------------------------------------
% METRICS
% -------------------------------------------------------------------------

% Get coverage of gray-matter a priori mask space (do we have data in all vox?)
% -------------------------------------------------------------------------

objmasked = apply_mask(obj, mean(gwcsf{1}));

wh = isnan(objmasked.dat) | objmasked.dat == 0;
individual_metrics.gray_matter_coverage = sum(~wh) ./ size(objmasked.dat, 1);
group_metrics.mean_gray_matter_coverage = mean(individual_metrics.gray_matter_coverage);


%brightedge
% NOT DONE YET

% Global activation/deactivation in CSF and white matter
% -------------------------------------------------------------------------

individual_metrics.global_ventricle_values = mean(values(:, 3));

group_metrics.global_d_ventricles = mean(values(:, 3)) ./ std(values(:, 3));

n = size(values, 1);
t = group_metrics.global_d_ventricles .* sqrt(n);   % t-value
p = 2 .* min(tcdf(t, n - 1), tcdf(-t, n - 1));      % two-tailed

group_metrics.global_logp_ventricles = -log(p) ./ (-log(.05));

% --- wm ---
individual_metrics.global_white_matter_values = mean(values(:, 2));

group_metrics.global_d_wm = mean(values(:, 2)) ./ std(values(:, 2));

t = group_metrics.global_d_wm .* sqrt(n);   % t-value
p = 2 .* min(tcdf(t, n - 1), tcdf(-t, n - 1));      % two-tailed

group_metrics.global_logp_wm = -log(p) ./ (-log(.05));


% Gray-matter signal explained by CSF
% -------------------------------------------------------------------------

[r, p] = corr(values(:, 1), values(:, 3));

group_metrics.gm_explained_by_csf_pvalue = p;
group_metrics.r2_explained_by_csf = r .^ 2;

% Gray-matter scale (l2 norm) explained by CSF scale
% -------------------------------------------------------------------------

[r, p] = corr(gwcsf_l2norm(:, 1), gwcsf_l2norm(:, 3));

group_metrics.gm_l2norm_explained_by_csf_pvalue = p;
group_metrics.r2_l2norm_explained_by_csf = r .^ 2;

% Signal in ventricles
% -------------------------------------------------------------------------

csfsignal = mean(abs(gwcsf{3}.dat)) ./ mean(abs(gwcsf{1}.dat));
group_metrics.csf_to_gm_signal_ratio = mean(csfsignal);
individual_metrics.csf_to_gm_signal_ratio = csfsignal;

% Scale Inhomogeneity
% -------------------------------------------------------------------------

L1norm = sum(abs(gwcsf{1}.dat))';
group_metrics.gm_scale_inhom = std(L1norm) ./ mean(L1norm);
individual_metrics.gm_L1norm = L1norm;

L1norm = sum(abs(gwcsf{3}.dat))';
group_metrics.csf_scale_inhom = std(L1norm) ./ mean(L1norm);
individual_metrics.csf_L1norm = L1norm;

% Add warnings
% -------------------------------------------------------------------------
group_metrics.warnings = {};

if any(individual_metrics.gray_matter_coverage < .9)
    group_metrics.warnings{end + 1} = sprintf('Warning: Incomplete coverage, missing gray-matter data in some images.\n   - Mean coverage is d = %3.2f%%\n', group_metrics.mean_gray_matter_coverage .* 100);
end

if group_metrics.global_logp_ventricles > 1
    group_metrics.warnings{end + 1} = sprintf('Warning: Significant global activation in CSF space/ventricles.\n   - Effect size is d = %3.2f\n', group_metrics.global_d_ventricles);
end

if group_metrics.global_logp_wm > 1
    group_metrics.warnings{end + 1} = sprintf('Warning: Significant global activation in white matter.\n   - Effect size is d = %3.2f\n', group_metrics.global_d_wm);
end

if group_metrics.gm_explained_by_csf_pvalue < .05
    group_metrics.warnings{end + 1} = sprintf('Warning: Gray-matter individual diffs significantly correlated with mean CSF value.\n   - Var explained (r^2) = %3.2f%%', group_metrics.r2_explained_by_csf .* 100);
end

if group_metrics.gm_l2norm_explained_by_csf_pvalue < .05
    group_metrics.warnings{end + 1} = sprintf('Warning: Gray-matter scale (L2 norm) significantly correlated with mean CSF L2 norm.\n   - Var explained (r^2) = %3.2f%%', group_metrics.r2_l2norm_explained_by_csf .* 100);
end

if group_metrics.csf_to_gm_signal_ratio > 1
    group_metrics.warnings{end + 1} = sprintf('Warning: Strong non-zero signal in CSF relative to gray matter.\n   - Ratio is = %3.2f', group_metrics.csf_to_gm_signal_ratio);
end

if group_metrics.csf_scale_inhom > .3
    group_metrics.warnings{end + 1} = sprintf('Warning: High individual diffs in image scaling estimated from CSF.\n   - CV is = %3.2f', group_metrics.csf_scale_inhom);
end



% Display warnings
% -------------------------------------------------------------------------

if doverbose
   
    disp(' ');
    disp(group_metrics)
    disp(char(group_metrics.warnings{:}))
    
end

end % function



