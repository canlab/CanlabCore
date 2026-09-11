function [results, maps] = canlab_effect_size_map(input_data, varargin)
% Voxel-wise effect sizes, within- vs. between-subject variance decomposition, and power maps for a group fMRI contrast
%
% This function takes a set of first-level (single-subject) contrast
% images, and optionally the matching first-level residual-variance
% (ResMS) images and design matrices, and produces:
%
%   1. Effect-size maps (Cohen's d) for the group contrast.
%   2. A decomposition of the variance of the contrast across subjects into
%      a within-subject (measurement noise, scaled by design efficiency)
%      component and a between-subject (true individual differences)
%      component.
%   3. Power maps for an exact replication of the study (same N, same
%      amount of data per subject), at uncorrected and family-wise-error
%      (FWE) corrected thresholds.
%   4. A "scan-time tradeoff" analysis: given a fixed budget of scanner
%      hours, how should you split it between more subjects and more scan
%      time per subject to maximize FWE-corrected power? The tradeoff
%      depends on the balance of within- and between-subject variance,
%      which is why the decomposition in (2) matters.
%
% The key model is the standard two-level summary-statistics model. For
% one subject, the contrast estimate is
%
%     cope_i = c' * beta_i,   Var(cope_i | subject) = sigma2_r_i * c'(X_i'X_i)^-1 c
%
% where sigma2_r_i is the first-level residual variance (the ResMS image)
% and c'(X'X)^-1 c is the "design variance" (design inefficiency) for the
% contrast. Across subjects, the observed variance of cope_i is
%
%     sigma2_total = sigma2_between + sigma2_within
%     sigma2_within = mean_i [ sigma2_r_i * c'(X_i'X_i)^-1 c ]
%
% so sigma2_between is estimated by subtraction (floored at 0). The
% within-subject component scales as 1/(number of images per subject), so
% for a proposed design with n_new images per subject we use
%
%     sigma2_within(n_new) = sigma2_within * n_scans / n_new
%
% and the expected group t-statistic (noncentrality parameter) is
%
%     t_expected = |cope| * sqrt(N) / sqrt(sigma2_between + sigma2_within(n_new))
%
% Power is then 1 - nctcdf(u, N - 1, t_expected) for a critical value u.
%
% :Usage:
% ::
%
%     % (a) From a folder of first-level SPM analyses, one subfolder per subject
%     [results, maps] = canlab_effect_size_map(spm_parent_dir, 'contrast', 2, [optional inputs])
%
%     % (b) From fmri_data objects (one image per subject in each)
%     [results, maps] = canlab_effect_size_map(con_obj, 'resms', resms_obj, 'X', X, 'c', c, [optional inputs])
%     [results, maps] = canlab_effect_size_map(con_obj, resms_obj, 'SPM', SPM, [optional inputs])
%     [results, maps] = canlab_effect_size_map(con_obj, [optional inputs])   % no within/between split
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2010, 2026 Tor Wager
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
%   **input_data:**
%        One of:
%
%        - A folder name (char or string). The folder is searched
%          recursively for SPM.mat files; each folder containing one is
%          treated as one subject's first-level analysis. The contrast
%          image (SPM.xCon(k).Vcon), residual-variance image (SPM.VResMS),
%          whitened/filtered design matrix, and contrast vector are read
%          from each. Requires the 'contrast' option.
%
%        - A cell array of subject folder names (each containing SPM.mat),
%          treated the same way.
%
%        - An fmri_data (or other image_vector) object whose images
%          (columns of .dat) are the first-level contrast images, one per
%          subject. Supply ResMS images and design information with the
%          optional inputs below. Images may also be given as a character
%          array or cell array of filenames.
%
% :Optional Inputs:
%
%   **resms_obj (2nd positional argument), or 'resms', [image_vector or filenames]:**
%        First-level residual-variance (ResMS) images, one per subject, in
%        the same subject order as the contrast images. Needed for the
%        within/between decomposition. Ignored in folder mode (they are
%        read from SPM.mat). If omitted, all variance is treated as
%        between-subject and the scan-time tradeoff becomes a plain
%        power-vs-N curve.
%
%   **'contrast', [integer or name]:**
%        Folder mode only. Index into SPM.xCon, or the contrast name
%        (SPM.xCon(k).name) to analyze. Required in folder mode.
%
%   **'SPM', [struct or path to SPM.mat]:**
%        Object mode. A representative first-level SPM structure used to
%        get the design matrix, contrast vector, number of scans, TR, and
%        effective df. Use with 'contrast' to pick the contrast.
%
%   **'X', [n_scans x k matrix] and 'c', [k x 1 vector]:**
%        Object mode. A representative first-level (whitened/filtered)
%        design matrix and contrast vector, used to compute the design
%        variance c'(X'X)^-1 c and the number of scans. If you have SPM.mat,
%        use SPM.xX.xKXs.X and SPM.xCon(k).c.
%
%   **'design_variance', [scalar or N x 1 vector]:**
%        Object mode. Supply c'(X'X)^-1 c directly, either one value for
%        all subjects or one per subject. Requires 'nscan'.
%
%   **'nscan', [scalar or N x 1 vector]:**
%        Number of image volumes per subject in the original first-level
%        design. Taken from SPM.nscan or size(X, 1) when available.
%
%   **'TR', [scalar]:**
%        Repetition time in seconds, for the scan-time tradeoff. Default:
%        SPM.xY.RT when available, otherwise 2.
%
%   **'mask', [filename or image_vector]:**
%        Restrict the analysis to this mask (e.g., a localizer, gray-matter
%        mask, or a-priori region). Default: all voxels with valid
%        (nonzero, finite) data in every subject.
%
%   **'group_con', [image] and 'group_t', [image]:**
%        Optional group-level contrast and t images (filenames or
%        image_vector objects), e.g., from a robust regression. If given,
%        the total across-subject variance is recovered as
%        cope^2 * N / t^2 instead of being computed directly from the
%        contrast images with OLS. Both must be supplied together.
%
%   **'hours', [scalar]:**
%        Total scanner hours available for the tradeoff analysis.
%        Default = 60.
%
%   **'hours_per_session', 'setup_min_first', 'setup_min_repeat', 'N_range', 'min_images', 'df_shrink_factor':**
%        Passed to canlab_scan_time_allocation; see its help. Defaults:
%        1.5 hours/session, 30 and 15 minutes setup, N_range = 3:100,
%        min_images = number of design columns (folder / X mode) or 30,
%        df_shrink_factor = mean(SPM.xX.erdf / nscan) when available,
%        otherwise 0.84.
%
%   **'alpha_corrected', [scalar]:**
%        FWE alpha for the corrected power maps. Default = 0.05.
%
%   **'alpha_uncorrected', [vector]:**
%        Uncorrected (two-tailed) alpha levels for the uncorrected power
%        curves. Default = [0.05 0.001].
%
%   **'correction', ['auto' | 'grf' | 'bonferroni']:**
%        How to obtain the FWE-corrected critical t value. 'grf' uses
%        random field theory (spm_uc_RF) with estimated smoothness;
%        'bonferroni' uses the number of in-mask voxels; 'auto' (default)
%        takes the smaller of the two thresholds, as SPM does, and falls
%        back to Bonferroni if SPM is not on the path.
%
%   **'fwhm', [1 x 3 vector, in voxels]:**
%        Skip smoothness estimation and use this FWHM (e.g., SPM.xVol.FWHM
%        from a group-level SPM.mat).
%
%   **'resels', [1 x 4 vector]:**
%        Skip smoothness estimation and use these resel counts (e.g.,
%        SPM.xVol.R). A scalar is treated as the 3-D resel count.
%
%   **'power_method', ['noncentral' | 'shift']:**
%        'noncentral' (default) computes power from the noncentral t
%        distribution (nctcdf), which is exact under the model. 'shift'
%        uses the faster classical approximation 1 - tcdf(u - t_expected, df).
%
%   **'outputdir', [folder name]:**
%        If given, write the maps as NIfTI images, a text log of the
%        printed report, and a .mat file with results to this folder.
%        Default: nothing is written to disk.
%
%   **'doplot', [logical flag]:**
%        Show orthviews of the variance/power maps and the allocation
%        figure. Default = true. 'noplot' to turn off.
%
%   **'verbose', [logical flag]:**
%        Print the report. Default = true. 'noverbose' to turn off.
%
% :Outputs:
%
%   **results:**
%        Structure with fields (voxel vectors are [n_voxels x 1], in-mask):
%
%        - .N, .n_voxels, .subject_names, .contrast_name
%        - .design         design variance, nscan, TR, columns, etc.
%        - .cope, .t, .d   group contrast, t, and Cohen's d per voxel
%        - .sig2_total, .sig2_residual, .sig2_within, .sig2_between
%                          variance components per voxel
%        - .unit_sig2_within  within component for one image per subject
%        - .summary        means and percentiles of the above
%        - .smoothness     FWHM, resel counts, NISC (effective number of
%                          independent spatial comparisons), method
%        - .replication    power for an exact replication (same N and
%                          data per subject): .u_corrected,
%                          .power_corrected, .power_uncorrected
%        - .allocation     scan-time tradeoff: candidate N and per-N
%                          session/time/image counts (from
%                          canlab_scan_time_allocation), .expected_t
%                          [voxels x N], .power_corrected [voxels x N],
%                          .power_uncorrected [voxels x N x alphas],
%                          .u_corrected, .mean_power_*, and .best with the
%                          allocation that maximizes mean corrected power
%        - .power_map      corrected power per voxel at the best allocation
%        - .options        the option values used
%
%   **maps:**
%        Structure of fmri_data objects for the key maps, ready for
%        orthviews / montage / write: .cohens_d, .residual_std,
%        .within_subjects_std, .between_subjects_std, .power_replication,
%        .expected_t_best, .power_best.
%
% :Examples:
% ::
%
%    % ---------------------------------------------------------------
%    % Example 1: Synthetic data, no files required
%    % ---------------------------------------------------------------
%    % Use the 30 sample contrast images as "first-level" cope images and
%    % make up ResMS images and a first-level design to go with them.
%    con_obj = load_image_set('emotionreg', 'noverbose');
%    N = size(con_obj.dat, 2);
%    resms_obj = con_obj;
%    resms_obj.dat = 20 + 5 * rand(size(con_obj.dat));  % positive residual variances
%    X = [randn(200, 1) ones(200, 1)];                  % 200 scans, 2 columns
%    c = [1 0]';
%    [results, maps] = canlab_effect_size_map(con_obj, resms_obj, 'X', X, 'c', c, ...
%        'hours', 60, 'TR', 2, 'correction', 'bonferroni', 'noplot');
%    orthviews(maps.power_best);
%
%    % ---------------------------------------------------------------
%    % Example 2: A folder of first-level SPM analyses
%    % ---------------------------------------------------------------
%    % spm_dir/sub-01/SPM.mat, spm_dir/sub-02/SPM.mat, ... Contrast #2 is
%    % the one to analyze. Restrict to a gray-matter mask and write maps.
%    [results, maps] = canlab_effect_size_map('/data/study/first_level', ...
%        'contrast', 2, 'mask', which('gray_matter_mask.img'), ...
%        'hours', 100, 'outputdir', '/data/study/effect_size_maps');
%
%    % Report the best split of scanner hours
%    disp(results.allocation.best)
%
%    % ---------------------------------------------------------------
%    % Example 3: Robust-regression group maps + subject-level images
%    % ---------------------------------------------------------------
%    [results, maps] = canlab_effect_size_map(con_obj, resms_obj, 'SPM', SPM, ...
%        'contrast', 2, 'group_con', 'rob_beta_0001.img', 'group_t', 'rob_tmap_0001.img');
%
% :References:
%   Mumford, J. A. & Nichols, T. E. (2008). Power calculation for group
%   fMRI studies accounting for arbitrary design and temporal
%   autocorrelation. NeuroImage, 39, 261-268.
%
%   Wager, T. D., Lindquist, M. A., et al. (2009). Brain mediators of
%   cardiovascular responses to social threat, Part I. NeuroImage, 47,
%   821-835.
%
% :See also:
%   - canlab_scan_time_allocation, canlab_power_allocation_curves,
%     power_from_variance, fmri_data.ttest, statistic_image.threshold
%

% ..
%    Programmers' notes:
%    2010: Tor Wager. Original effect_size_map.m script, used for power
%          examples in Wager & Lindquist teaching materials.
%    2026-09: Tor Wager. Converted to a function with inputParser; accepts
%          SPM folders or fmri_data objects. Fixes relative to the script:
%          - The loop computing the p-value thresholds for NISC overwrote
%            the same variable on every iteration (only the last N's df
%            was used); it now stores one value per N.
%          - Uncorrected critical t values used df = N; now df = N - 1,
%            matching the one-sample t-test used for power.
%          - Smoothness was estimated from first-level ResMS images, which
%            are not residuals. It is now estimated from the standardized
%            residuals of the group one-sample t-test (cope images minus
%            their mean), which is what random field theory needs.
%          - RPV.img was read back from disk after spm_est_smoothness; SPM12+
%            writes RPV.nii. The resel counts are now taken from the third
%            output of spm_est_smoothness instead.
%          - The df shrink factor was applied to proposed designs but not
%            to the observed design when rescaling within-subject variance,
%            inflating within-subject variance for proposed designs by
%            ~19%. Scaling now uses raw scan counts on both sides; the
%            shrink factor is only used for the minimum-data feasibility check.
%          - Power used the shifted-central-t approximation; the default is
%            now the noncentral t distribution ('shift' is still available).
%          - Per-subject design variance and scan counts are used when
%            available (folder mode) instead of one example SPM.mat.
%          - No more cd() into an output folder or unconditional diary().
% ..

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

% Positional ResMS object: canlab_effect_size_map(con_obj, resms_obj, ...)
resms_in = [];
if ~isempty(varargin) && (isa(varargin{1}, 'image_vector') || (isnumeric(varargin{1}) && isempty(varargin{1})))
    resms_in = varargin{1};
    varargin(1) = [];
end

% Parse special command keywords and remove them before inputParser

doplot = true;
plot_idx = strcmpi(varargin, 'noplot');
if any(plot_idx)
    doplot = false;
    varargin(plot_idx) = [];   % remove so inputParser doesn't see it
end
plot_idx = strcmpi(varargin, 'plot');
if any(plot_idx)               % Override: omit 'doplot' key/value pair
    doplot = true;
    varargin(plot_idx) = [];
end

verbose = true;
verbose_idx = strcmpi(varargin, 'noverbose');
if any(verbose_idx)
    verbose = false;
    varargin(verbose_idx) = [];   % remove so inputParser doesn't see it
end

% Use inputParser to parse key/value pairs
% First add obligatory/non-conditional keywords

valfcn_image = @(x) isempty(x) || isa(x, 'image_vector') || ischar(x) || isstring(x) || iscell(x);
valfcn_logical = @(x) islogical(x) || isnumeric(x);
valfcn_posscalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', 'positive'});
valfcn_scalar_or_empty = @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0);

p = inputParser;
p.PartialMatching = false;   % 'c' and 'contrast' are both parameter names
p.addRequired('input_data', @(x) isa(x, 'image_vector') || ischar(x) || isstring(x) || iscell(x));
p.addParameter('resms', resms_in, valfcn_image);
p.addParameter('contrast', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)) || ischar(x) || isstring(x));
p.addParameter('SPM', [], @(x) isempty(x) || isstruct(x) || ischar(x) || isstring(x));
p.addParameter('X', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
p.addParameter('c', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('design_variance', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('nscan', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('TR', [], valfcn_scalar_or_empty);
p.addParameter('mask', [], valfcn_image);
p.addParameter('group_con', [], valfcn_image);
p.addParameter('group_t', [], valfcn_image);

p.addParameter('hours', 60, valfcn_posscalar);
p.addParameter('hours_per_session', 1.5, valfcn_posscalar);
p.addParameter('setup_min_first', 30, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
p.addParameter('setup_min_repeat', 15, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
p.addParameter('N_range', 3:100, @(x) validateattributes(x, {'numeric'}, {'nonempty', 'vector', 'positive', 'integer'}));
p.addParameter('min_images', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
p.addParameter('df_shrink_factor', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0 && x <= 1));

p.addParameter('alpha_corrected', 0.05, @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0, '<', 1}));
p.addParameter('alpha_uncorrected', [0.05 0.001], @(x) validateattributes(x, {'numeric'}, {'nonempty', 'vector', '>', 0, '<', 1}));
p.addParameter('correction', 'auto', @(x) any(strcmpi(x, {'auto', 'grf', 'bonferroni'})));
p.addParameter('fwhm', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 3));
p.addParameter('resels', [], @(x) isempty(x) || (isnumeric(x) && (numel(x) == 1 || numel(x) == 4)));
p.addParameter('power_method', 'noncentral', @(x) any(strcmpi(x, {'noncentral', 'shift'})));
p.addParameter('outputdir', '', @(x) ischar(x) || isstring(x));

% Special key/value pairs that we have potentially set with optional keywords
p.addParameter('doplot', doplot, valfcn_logical);
p.addParameter('verbose', verbose, valfcn_logical);

% process inputs and deal out to variables in workspace
p.parse(input_data, varargin{:});
ARGS = p.Results;

doplot = logical(ARGS.doplot);
verbose = logical(ARGS.verbose);
outputdir = char(ARGS.outputdir);
correction = lower(char(ARGS.correction));
power_method = lower(char(ARGS.power_method));
alpha_c = ARGS.alpha_corrected;
alpha_unc = ARGS.alpha_uncorrected(:)';

if xor(isempty(ARGS.group_con), isempty(ARGS.group_t))
    error('canlab_effect_size_map:GroupMaps', '''group_con'' and ''group_t'' must be supplied together.');
end

% Display helpers
dashes = '------------------------------------------------------------';
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

% -------------------------------------------------------------------------
% Output folder and log file
% -------------------------------------------------------------------------

if ~isempty(outputdir)
    if ~exist(outputdir, 'dir'), mkdir(outputdir); end
    logfile = fullfile(outputdir, 'canlab_effect_size_map_log.txt');
    if verbose, fprintf('Saving text output in %s\n', logfile); end
    diary(logfile);
    diary_cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
end

if verbose, printhdr('canlab_effect_size_map: loading data'); end

% -------------------------------------------------------------------------
% Load data: folder mode or object mode
% -------------------------------------------------------------------------

design = struct('desvar', [], 'nscan', [], 'k', [], 'erdf', [], 'TR', [], 'contrast_name', '', 'c', [], 'source', '');

is_folder_mode = (ischar(input_data) || isstring(input_data) || iscell(input_data)) && looks_like_folders(input_data);

if is_folder_mode
    % ---- Folder mode: one SPM.mat per subject ---------------------------
    if isempty(ARGS.contrast)
        error('canlab_effect_size_map:NoContrast', 'Folder mode requires the ''contrast'' option (index into SPM.xCon or contrast name).');
    end

    [con_obj, resms_obj, design, subject_names] = load_from_spm_folders(input_data, ARGS.contrast, verbose);

else
    % ---- Object mode -----------------------------------------------------
    con_obj = to_fmri_data(input_data, 'contrast images');
    subject_names = cellstr(con_obj.image_names);

    resms_obj = [];
    if ~isempty(ARGS.resms)
        resms_obj = to_fmri_data(ARGS.resms, 'ResMS images');
    end

    design = design_from_options(ARGS, size(con_obj.dat, 2));
end

N = size(con_obj.dat, 2);

if N < 3
    error('canlab_effect_size_map:TooFewSubjects', 'At least 3 subjects (contrast images) are needed; found %d.', N);
end

if ~isempty(resms_obj) && size(resms_obj.dat, 2) ~= N
    error('canlab_effect_size_map:ResMSMismatch', 'Number of ResMS images (%d) does not match number of contrast images (%d).', size(resms_obj.dat, 2), N);
end

% Fill in TR / df shrink / min_images defaults from the design if not given
if ~isempty(ARGS.TR), design.TR = ARGS.TR; end
if isempty(design.TR), design.TR = 2; end

df_shrink_factor = ARGS.df_shrink_factor;
if isempty(df_shrink_factor)
    if ~isempty(design.erdf) && ~isempty(design.nscan)
        df_shrink_factor = mean(design.erdf(:) ./ design.nscan(:), 'omitnan');
        df_shrink_factor = min(max(df_shrink_factor, eps), 1);
    else
        df_shrink_factor = 0.84;
    end
end

min_images = ARGS.min_images;
if isempty(min_images)
    if ~isempty(design.k), min_images = design.k; else, min_images = 30; end
end

% -------------------------------------------------------------------------
% Mask and valid voxels
% -------------------------------------------------------------------------

if ~isempty(ARGS.mask)
    if verbose, fprintf('Applying mask\n'); end
    con_obj = apply_mask(con_obj, ARGS.mask);
    if ~isempty(resms_obj), resms_obj = apply_mask(resms_obj, ARGS.mask); end
end

% Put everything in the same voxel space as the contrast images
con_obj = replace_empty(con_obj);

if ~isempty(resms_obj)
    if compare_space(resms_obj, con_obj)
        if verbose, fprintf('Resampling ResMS images to contrast image space\n'); end
        resms_obj = resample_space(resms_obj, con_obj);
    end
    resms_obj = replace_empty(resms_obj);
    if size(resms_obj.dat, 1) ~= size(con_obj.dat, 1)
        error('canlab_effect_size_map:SpaceMismatch', 'ResMS and contrast images do not have the same number of in-mask voxels after resampling.');
    end
end

[gcon_in, gt_in] = deal([]);
if ~isempty(ARGS.group_con)
    gcon_in = load_in_space(ARGS.group_con, con_obj, 'group_con');
    gt_in = load_in_space(ARGS.group_t, con_obj, 'group_t');
end

% Voxels must have valid data in every subject (and in the group maps)
wh_bad = any(~isfinite(con_obj.dat) | con_obj.dat == 0, 2);
if ~isempty(resms_obj)
    wh_bad = wh_bad | any(~isfinite(resms_obj.dat) | resms_obj.dat <= 0, 2);
end
if ~isempty(gcon_in)
    wh_bad = wh_bad | ~isfinite(gcon_in) | ~isfinite(gt_in) | gt_in == 0;
end

con_obj = remove_empty(con_obj, wh_bad);
if ~isempty(resms_obj), resms_obj = remove_empty(resms_obj, wh_bad); end
if ~isempty(gcon_in), gcon_in = gcon_in(~wh_bad); gt_in = gt_in(~wh_bad); end

nvox = size(con_obj.dat, 1);
if nvox == 0
    error('canlab_effect_size_map:NoVoxels', 'No voxels with valid data in all subjects.');
end
if size(con_obj.dat, 2) ~= N || (~isempty(resms_obj) && size(resms_obj.dat, 2) ~= N)
    error('canlab_effect_size_map:EmptyImages', 'One or more contrast or ResMS images contain no valid data (all zero or NaN).');
end

% A one-image template object we can copy to create output maps
template = get_wh_image(con_obj, 1);
template.image_names = '';
template.fullpath = '';
template.history = {};
template.X = []; template.Y = [];
if isprop(template, 'metadata_table'), template.metadata_table = table(); end

if verbose
    fprintf('Subjects: N = %d\n', N);
    fprintf('Voxels with valid data in all subjects: %d\n', nvox);
    if ~isempty(design.contrast_name), fprintf('Contrast: %s\n', design.contrast_name); end
end

% -------------------------------------------------------------------------
% Group statistics and variance components (all [nvox x 1])
% -------------------------------------------------------------------------

con = double(con_obj.dat);                       % voxels x N

if ~isempty(gcon_in)
    % Group maps supplied (e.g., robust regression): recover the implied
    % total variance from the t-statistic. t = cope * sqrt(N) / sigma.
    cope = double(gcon_in(:));
    t = double(gt_in(:));
    sig2_total = cope .^ 2 .* N ./ t .^ 2;
    group_source = 'supplied group_con / group_t images';
else
    % OLS one-sample t-test across subjects
    cope = mean(con, 2);
    sd = std(con, 0, 2);
    t = cope ./ (sd ./ sqrt(N));
    sig2_total = sd .^ 2;
    group_source = 'one-sample t-test on contrast images';
end

d = t ./ sqrt(N);                                 % Cohen's d (signed)

% Within-subject component: residual variance * design variance, per
% subject, then averaged. unit_sig2w is the contribution of a single image
% volume, i.e., sig2w * nscan, which lets us rescale to other scan counts.
have_within = ~isempty(resms_obj) && ~isempty(design.desvar);

if have_within
    resms = double(resms_obj.dat);                % voxels x N
    desvar = expand_per_subject(design.desvar, N, 'design_variance');
    nscan = expand_per_subject(design.nscan, N, 'nscan');

    sig2_residual = mean(resms, 2);
    sig2_within = mean(resms .* desvar', 2);
    unit_sig2_within = mean(resms .* (desvar .* nscan)', 2);
else
    if ~isempty(resms_obj) && isempty(design.desvar) && verbose
        fprintf('ResMS images given but no design information (''SPM'', ''X''/''c'', or ''design_variance''): skipping within/between decomposition.\n');
    end
    [sig2_residual, sig2_within, unit_sig2_within] = deal(zeros(nvox, 1));
    desvar = []; nscan = design.nscan;
end

sig2_between = sig2_total - sig2_within;
n_floored = sum(sig2_between < 0);
sig2_between(sig2_between < 0) = 0;              % variances cannot be negative

% -------------------------------------------------------------------------
% Smoothness / resels, for FWE-corrected thresholds
% -------------------------------------------------------------------------

smoothness = get_smoothness(con, template, ARGS.fwhm, ARGS.resels, correction, verbose);

% -------------------------------------------------------------------------
% Power for an exact replication (same N, same data per subject)
% -------------------------------------------------------------------------

df_obs = N - 1;
u_corr_obs = corrected_threshold(alpha_c, df_obs, nvox, smoothness);
ncp_obs = abs(cope) .* sqrt(N) ./ sqrt(sig2_total);   % expected |t| under the fitted model
ncp_obs(~isfinite(ncp_obs)) = NaN;

replication = struct();
replication.N = N;
replication.expected_t = ncp_obs;
replication.u_corrected = u_corr_obs;
replication.power_corrected = power_from_ncp(ncp_obs, df_obs, u_corr_obs, power_method);
replication.power_uncorrected = zeros(nvox, numel(alpha_unc));
for a = 1:numel(alpha_unc)
    u_unc = tinv(1 - alpha_unc(a) / 2, df_obs);   % two-tailed
    replication.power_uncorrected(:, a) = power_from_ncp(ncp_obs, df_obs, u_unc, power_method);
end
replication.alpha_uncorrected = alpha_unc;

% -------------------------------------------------------------------------
% Scan-time tradeoff: candidate N given fixed scanner hours
% -------------------------------------------------------------------------

alloc = canlab_scan_time_allocation(ARGS.hours, 'N_range', ARGS.N_range, ...
    'hours_per_session', ARGS.hours_per_session, 'TR', design.TR, ...
    'setup_min_first', ARGS.setup_min_first, 'setup_min_repeat', ARGS.setup_min_repeat, ...
    'df_shrink_factor', df_shrink_factor, 'min_images', min_images);

nN = numel(alloc.N);
[expected_t, power_corrected] = deal(zeros(nvox, nN));
power_uncorrected = zeros(nvox, nN, numel(alpha_unc));
u_corrected = zeros(1, nN);

for i = 1:nN
    Ni = alloc.N(i);
    dfi = Ni - 1;

    % Within-subject variance for the proposed number of images per subject
    sig2wi = unit_sig2_within ./ alloc.images_per_subject(i);

    % Expected t-value (noncentrality parameter), unsigned
    ncp = abs(cope) .* sqrt(Ni) ./ sqrt(sig2_between + sig2wi);
    ncp(~isfinite(ncp)) = NaN;
    expected_t(:, i) = ncp;

    u_corrected(i) = corrected_threshold(alpha_c, dfi, nvox, smoothness);
    power_corrected(:, i) = power_from_ncp(ncp, dfi, u_corrected(i), power_method);

    for a = 1:numel(alpha_unc)
        u_unc = tinv(1 - alpha_unc(a) / 2, dfi);
        power_uncorrected(:, i, a) = power_from_ncp(ncp, dfi, u_unc, power_method);
    end
end

% Effective number of independent spatial comparisons (NISC): the number
% of independent tests you would Bonferroni-correct for to get the same
% corrected threshold. One value per N; averaged over N > 20, where the
% df-dependence of the threshold has settled down.
p_thresh = tcdf(u_corrected, alloc.N - 1, 'upper');
NISC_by_N = alpha_c ./ p_thresh;
wh_stable = alloc.N > 20;
if ~any(wh_stable), wh_stable = true(size(alloc.N)); end
smoothness.NISC = mean(NISC_by_N(wh_stable));
smoothness.NISC_by_N = NISC_by_N;

% Best allocation: maximize mean corrected power across the search area
mean_power_corrected = mean(power_corrected, 1, 'omitnan');
mean_power_uncorrected = reshape(mean(power_uncorrected, 1, 'omitnan'), nN, numel(alpha_unc));

[~, wh_best] = max(mean_power_corrected);

best = struct();
best.N = alloc.N(wh_best);
best.sessions_per_subject = alloc.sessions_per_subject(wh_best);
best.hours_per_subject = alloc.hours_per_subject(wh_best);
best.functional_hours_per_subject = alloc.functional_hours_per_subject(wh_best);
best.functional_min_per_subject = alloc.functional_min_per_subject(wh_best);
best.images_per_subject = alloc.images_per_subject(wh_best);
best.mean_power_corrected = mean_power_corrected(wh_best);
best.sig2_within_at_best = unit_sig2_within ./ alloc.images_per_subject(wh_best);
best.index = wh_best;

power_map = power_corrected(:, wh_best);

% -------------------------------------------------------------------------
% Summary statistics
% -------------------------------------------------------------------------

pct = [25 50 75];
summary_stats = struct();
summary_stats.mean_d = mean(d, 'omitnan');
summary_stats.mean_abs_d = mean(abs(d), 'omitnan');
summary_stats.prctile_d = prctile(d, pct);
summary_stats.mean_sig2_residual = mean(sig2_residual, 'omitnan');
summary_stats.mean_sig2_within = mean(sig2_within, 'omitnan');
summary_stats.mean_sig2_between = mean(sig2_between, 'omitnan');
summary_stats.mean_sig2_total = mean(sig2_total, 'omitnan');
summary_stats.prctile_sig2_within = prctile(sig2_within, pct);
summary_stats.prctile_sig2_between = prctile(sig2_between, pct);
summary_stats.mean_within_std_per_image = mean(sqrt(unit_sig2_within), 'omitnan');
summary_stats.mean_between_std = mean(sqrt(sig2_between), 'omitnan');
summary_stats.n_voxels_between_floored = n_floored;
summary_stats.percentiles = pct;

% -------------------------------------------------------------------------
% Report
% -------------------------------------------------------------------------

if verbose
    printhdr('Effect size and variance components');
    fprintf('Group statistics from: %s\n', group_source);
    fprintf('Mask area: %d voxels\n', nvox);
    fprintf('Subjects: N = %d\n', N);
    if ~isempty(design.nscan)
        fprintf('Scans per subject (first level): mean %3.1f\n', mean(design.nscan));
    end
    if ~isempty(desvar)
        fprintf('Design variance c''(X''X)^-1 c: mean %3.5f (per-scan: %3.4f)\n', mean(desvar), mean(desvar .* nscan));
    end
    fprintf('df shrink factor (effective df / scans): %3.3f\n', df_shrink_factor);
    fprintf('\nPooled average effect size (Cohen''s d): %3.2f  (|d|: %3.2f)\n', summary_stats.mean_d, summary_stats.mean_abs_d);
    fprintf('Cohen''s d percentiles (25/50/75): %3.2f / %3.2f / %3.2f\n', summary_stats.prctile_d);

    if have_within
        fprintf('\nVariance components (mean across voxels):\n');
        fprintf('  Residual variance, first level (sig2_residual):         %3.3f\n', summary_stats.mean_sig2_residual);
        fprintf('  Within-subject contribution to cope variance (sig2_w):  %3.3f\n', summary_stats.mean_sig2_within);
        fprintf('  Between-subject variance (sig2_b):                      %3.3f\n', summary_stats.mean_sig2_between);
        fprintf('  Total cope variance across subjects (sig2_total):       %3.3f\n', summary_stats.mean_sig2_total);
        fprintf('  Voxels where sig2_b was floored at 0: %d (%3.1f%%)\n', n_floored, 100 * n_floored / nvox);
        fprintf('\nStandard deviations, so that se(cope) = sqrt(std_w^2 / nscan + std_b^2 / N):\n');
        fprintf('  Within-subject std per image volume (std_w): %3.3f\n', summary_stats.mean_within_std_per_image);
        fprintf('  Between-subject std (std_b):                 %3.3f\n', summary_stats.mean_between_std);
        fprintf('\nVariance summary for in-mask voxels:\n');
        fprintf('Percentile:   \t25th\t50th\t75th\n');
        fprintf('Within-subj:  \t%3.2f\t%3.2f\t%3.2f\n', summary_stats.prctile_sig2_within);
        fprintf('Between-subj: \t%3.2f\t%3.2f\t%3.2f\n', summary_stats.prctile_sig2_between);
    else
        fprintf('\nNo within-subject variance information: all variance treated as between-subject.\n');
        fprintf('Total cope variance across subjects (mean): %3.3f\n', summary_stats.mean_sig2_total);
    end

    printhdr('Multiple comparisons');
    fprintf('Correction: %s\n', smoothness.description);
    if ~isempty(smoothness.fwhm)
        fprintf('Smoothness (FWHM, voxels): %3.1f %3.1f %3.1f\n', smoothness.fwhm);
    end
    if ~isempty(smoothness.resels)
        fprintf('Resel counts: %s\n', num2str(smoothness.resels, '%3.1f '));
    end
    fprintf('Effective number of independent spatial comparisons (NISC): %3.0f\n', smoothness.NISC);

    printhdr(sprintf('Power for an exact replication (N = %d)', N));
    fprintf('Critical t (FWE %3.3f, one-tailed): %3.2f\n', alpha_c, u_corr_obs);
    fprintf('FWE-corrected power: mean %3.0f%%, min %3.0f%%, max %3.0f%%, voxels >= 80%%: %d\n', ...
        100 * mean(replication.power_corrected, 'omitnan'), 100 * min(replication.power_corrected), ...
        100 * max(replication.power_corrected), sum(replication.power_corrected >= .8));
    for a = 1:numel(alpha_unc)
        fprintf('Uncorrected power at p < %g (two-tailed): mean %3.0f%%\n', alpha_unc(a), 100 * mean(replication.power_uncorrected(:, a), 'omitnan'));
    end

    printhdr(sprintf('Best balance of subjects and time for %3.0f total scanner hours', ARGS.hours));
    fprintf('Subjects: %d\n', best.N);
    fprintf('Sessions per subject (max session time of %3.1f hours): %d\n', ARGS.hours_per_session, best.sessions_per_subject);
    fprintf('Total scan hours per subject: %3.1f; functional time: %3.1f hours, or %3.0f min (%3.0f images at TR = %3.2f s)\n', ...
        best.hours_per_subject, best.functional_hours_per_subject, best.functional_min_per_subject, best.images_per_subject, design.TR);
    fprintf('\nPower for FWE-corrected search at the best allocation:\n');
    fprintf('Minimum within search area: %3.0f%%\n', 100 * min(power_map));
    fprintf('Maximum within search area: %3.0f%%\n', 100 * max(power_map));
    fprintf('Mean within search area:    %3.0f%%\n', 100 * mean(power_map, 'omitnan'));
    fprintf('Voxels with 80%% power:      %d\n', sum(power_map >= .8));
    if ~isempty(alloc.infeasible_N)
        fprintf('\n(N >= %d dropped: too little functional time per subject for a first-level model.)\n', min(alloc.infeasible_N));
    end
end

% -------------------------------------------------------------------------
% Assemble outputs
% -------------------------------------------------------------------------

results = struct();
results.N = N;
results.n_voxels = nvox;
results.subject_names = subject_names;
results.contrast_name = design.contrast_name;
results.group_source = group_source;

design.desvar_per_subject = desvar;
design.nscan_per_subject = nscan;
design.df_shrink_factor = df_shrink_factor;
design.min_images = min_images;
results.design = design;

results.cope = cope;
results.t = t;
results.d = d;
results.sig2_total = sig2_total;
results.sig2_residual = sig2_residual;
results.sig2_within = sig2_within;
results.sig2_between = sig2_between;
results.unit_sig2_within = unit_sig2_within;
results.have_within_decomposition = have_within;
results.summary = summary_stats;
results.smoothness = smoothness;
results.replication = replication;

allocation = alloc;
allocation.expected_t = expected_t;
allocation.power_corrected = power_corrected;
allocation.power_uncorrected = power_uncorrected;
allocation.alpha_corrected = alpha_c;
allocation.alpha_uncorrected = alpha_unc;
allocation.u_corrected = u_corrected;
allocation.mean_power_corrected = mean_power_corrected;
allocation.mean_power_uncorrected = mean_power_uncorrected;
allocation.best = best;
results.allocation = allocation;

results.power_map = power_map;
results.volInfo = con_obj.volInfo;
results.removed_voxels = con_obj.removed_voxels;
% Options used, without the (possibly large) data inputs
results.options = rmfield(ARGS, {'input_data', 'resms', 'group_con', 'group_t', 'mask', 'SPM', 'X'});

% Map objects
maps = struct();
maps.cohens_d = make_map(template, d, 'Cohen''s d');
maps.residual_std = make_map(template, sqrt(sig2_residual), 'First-level residual std');
maps.within_subjects_std = make_map(template, sqrt(sig2_within), 'Within-subject std (contribution to cope)');
maps.between_subjects_std = make_map(template, sqrt(sig2_between), 'Between-subject std');
maps.power_replication = make_map(template, replication.power_corrected, sprintf('FWE-corrected power, replication N = %d', N));
maps.expected_t_best = make_map(template, expected_t(:, wh_best), sprintf('Expected t, N = %d', best.N));
maps.power_best = make_map(template, power_map, sprintf('FWE-corrected power, N = %d, %d min/subject', best.N, best.functional_min_per_subject));

% -------------------------------------------------------------------------
% Write images and results
% -------------------------------------------------------------------------

if ~isempty(outputdir)
    pmapname = sprintf('power_map_%03dhours_N%03d_ftime%03dmin.nii', round(ARGS.hours), best.N, best.functional_min_per_subject);
    fnames = {'cohens_d.nii', 'residual_std.nii', 'within_subjects_std.nii', 'between_subjects_std.nii', ...
        sprintf('power_replication_N%03d.nii', N), 'expected_t_best.nii', pmapname};
    mapnames = fieldnames(maps);

    for i = 1:numel(mapnames)
        write_map(maps.(mapnames{i}), fullfile(outputdir, fnames{i}), verbose);
    end

    save(fullfile(outputdir, 'canlab_effect_size_map_results.mat'), 'results', 'maps');
    if verbose, fprintf('Saved results to %s\n', fullfile(outputdir, 'canlab_effect_size_map_results.mat')); end
end

% -------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------

if doplot
    plot_maps(maps, have_within);
    plot_allocation(results, have_within);

    if ~isempty(outputdir)
        try
            scn_export_papersetup(600);
            saveas(gcf, fullfile(outputdir, 'Optimal_allocation_of_hours'), 'png');
            saveas(gcf, fullfile(outputdir, 'Optimal_allocation_of_hours'), 'fig');
        catch
            disp('Cannot save allocation figure')
        end
    end
end

end % main function



% =========================================================================
% Subfunctions: input handling
% =========================================================================

function tf = looks_like_folders(x)
% True if x names one or more folders (SPM folder mode) rather than image files.
if iscell(x)
    tf = all(cellfun(@(s) (ischar(s) || isstring(s)) && exist(char(s), 'dir') == 7, x));
else
    x = char(x);
    tf = size(x, 1) == 1 && exist(x, 'dir') == 7;
end
end


function obj = to_fmri_data(x, what)
% Convert filenames or image_vector objects to fmri_data
if isa(x, 'fmri_data')
    obj = x;
elseif isa(x, 'image_vector')
    obj = fmri_data(x);
elseif ischar(x) || isstring(x) || iscell(x)
    x = char(x);
    obj = fmri_data(x, 'noverbose');
else
    error('canlab_effect_size_map:BadInput', 'Cannot interpret %s input.', what);
end
end


function vals = load_in_space(x, con_obj, what)
% Load a group map and return its voxel values in con_obj's (full) space
obj = to_fmri_data(x, what);
if size(obj.dat, 2) ~= 1
    error('canlab_effect_size_map:GroupMaps', '%s must be a single image.', what);
end
if compare_space(obj, con_obj)
    obj = resample_space(obj, con_obj);
end
obj = replace_empty(obj);
if size(obj.dat, 1) ~= size(con_obj.dat, 1)
    error('canlab_effect_size_map:SpaceMismatch', '%s does not match the contrast image space after resampling.', what);
end
vals = double(obj.dat);
end


function v = expand_per_subject(v, N, what)
% Return an N x 1 vector from a scalar or N-vector
v = double(v(:));
if isscalar(v)
    v = repmat(v, N, 1);
elseif numel(v) ~= N
    error('canlab_effect_size_map:BadLength', '''%s'' must be a scalar or have one value per subject (%d).', what, N);
end
end


function design = design_from_options(ARGS, N)
% Build the design struct in object mode from 'SPM', 'X'/'c', or 'design_variance'
design = struct('desvar', [], 'nscan', [], 'k', [], 'erdf', [], 'TR', [], 'contrast_name', '', 'c', [], 'source', '');

if ~isempty(ARGS.SPM)
    SPM = ARGS.SPM;
    if ischar(SPM) || isstring(SPM)
        S = load(char(SPM), 'SPM');
        SPM = S.SPM;
    end
    design = design_from_spm(SPM, ARGS.contrast);
    design.source = 'SPM';

elseif ~isempty(ARGS.X)
    if isempty(ARGS.c)
        error('canlab_effect_size_map:NoContrastVector', '''X'' requires a contrast vector ''c''.');
    end
    X = double(ARGS.X);
    c = double(ARGS.c(:));
    if numel(c) ~= size(X, 2)
        error('canlab_effect_size_map:ContrastSize', '''c'' must have one element per column of ''X'' (%d).', size(X, 2));
    end
    pX = pinv(X);
    design.desvar = c' * (pX * pX') * c;
    design.nscan = size(X, 1);
    design.k = size(X, 2);
    design.c = c;
    design.source = 'X and c';

elseif ~isempty(ARGS.design_variance)
    design.desvar = ARGS.design_variance;
    if isempty(ARGS.nscan)
        error('canlab_effect_size_map:NoNscan', '''design_variance'' requires ''nscan'' (number of scans per subject).');
    end
    design.source = 'design_variance';
end

if ~isempty(ARGS.nscan), design.nscan = ARGS.nscan; end

if ~isempty(design.desvar)
    expand_per_subject(design.desvar, N, 'design_variance');  % validate length
end
end


function design = design_from_spm(SPM, contrast_spec)
% Extract design variance, scan counts, TR, and effective df from an SPM struct
design = struct('desvar', [], 'nscan', [], 'k', [], 'erdf', [], 'TR', [], 'contrast_name', '', 'c', [], 'source', 'SPM');

k = resolve_contrast(SPM, contrast_spec);
c = double(SPM.xCon(k).c);
if size(c, 2) > 1
    warning('canlab_effect_size_map:FContrast', 'Contrast %d has %d columns (F-contrast?). Using the first column.', k, size(c, 2));
    c = c(:, 1);
end

if isfield(SPM.xX, 'Bcov') && ~isempty(SPM.xX.Bcov)
    % Covariance of parameter estimates (unscaled), as used by SPM's spm_contrasts
    Bcov = SPM.xX.Bcov;
elseif isfield(SPM.xX, 'xKXs') && isfield(SPM.xX.xKXs, 'X')
    pX = pinv(SPM.xX.xKXs.X);
    Bcov = pX * pX';
else
    warning('canlab_effect_size_map:UnfilteredDesign', 'SPM.xX.Bcov / xKXs not found (model not estimated?). Using the raw design matrix SPM.xX.X.');
    pX = pinv(SPM.xX.X);
    Bcov = pX * pX';
end

design.desvar = c' * Bcov * c;
design.nscan = sum(SPM.nscan);
design.k = size(SPM.xX.X, 2);
design.c = c;
design.contrast_name = SPM.xCon(k).name;
design.contrast_index = k;

if isfield(SPM.xX, 'erdf') && ~isempty(SPM.xX.erdf), design.erdf = SPM.xX.erdf; end
if isfield(SPM, 'xY') && isfield(SPM.xY, 'RT') && ~isempty(SPM.xY.RT), design.TR = SPM.xY.RT; end
end


function k = resolve_contrast(SPM, contrast_spec)
% Contrast index from an index or a name
if ~isfield(SPM, 'xCon') || isempty(SPM.xCon)
    error('canlab_effect_size_map:NoContrasts', 'SPM.xCon is empty: no contrasts have been defined.');
end
if isempty(contrast_spec)
    error('canlab_effect_size_map:NoContrast', 'Specify which contrast to use with ''contrast'', [index or name].');
end
if isnumeric(contrast_spec)
    k = contrast_spec;
    if k < 1 || k > numel(SPM.xCon)
        error('canlab_effect_size_map:BadContrast', 'Contrast index %d is out of range (1-%d).', k, numel(SPM.xCon));
    end
else
    names = {SPM.xCon.name};
    k = find(strcmp(names, char(contrast_spec)));
    if isempty(k)
        error('canlab_effect_size_map:BadContrast', 'No contrast named ''%s''. Available: %s', char(contrast_spec), strjoin(names, ' | '));
    end
    k = k(1);
end
end


function [con_obj, resms_obj, design, subject_names] = load_from_spm_folders(input_data, contrast_spec, verbose)
% Find SPM.mat files under a parent folder (or in a list of folders) and
% load each subject's contrast image, ResMS image, and design information.

if iscell(input_data)
    subject_dirs = cellfun(@char, input_data(:), 'UniformOutput', false);
    wh_missing = ~cellfun(@(d) exist(fullfile(d, 'SPM.mat'), 'file') == 2, subject_dirs);
    if any(wh_missing)
        error('canlab_effect_size_map:NoSPM', 'No SPM.mat in: %s', strjoin(subject_dirs(wh_missing), ', '));
    end
else
    parent = char(input_data);
    dd = dir(fullfile(parent, '**', 'SPM.mat'));
    if isempty(dd)
        error('canlab_effect_size_map:NoSPM', 'No SPM.mat files found under %s', parent);
    end
    subject_dirs = unique({dd.folder}');
end

N = numel(subject_dirs);
if verbose, fprintf('Found %d first-level SPM analyses under %s\n', N, char(input_data(1, :))); end

[con_files, resms_files] = deal(cell(N, 1));
[desvar, nscan, erdf, TRs] = deal(nan(N, 1));
k_cols = nan(N, 1);
contrast_name = '';

for i = 1:N
    S = load(fullfile(subject_dirs{i}, 'SPM.mat'), 'SPM');
    SPM = S.SPM;

    dsn = design_from_spm(SPM, contrast_spec);
    k = dsn.contrast_index;

    if i == 1
        contrast_name = dsn.contrast_name;
    elseif ~strcmp(dsn.contrast_name, contrast_name)
        warning('canlab_effect_size_map:ContrastNameMismatch', ...
            'Contrast name differs across subjects: ''%s'' (%s) vs ''%s''.', dsn.contrast_name, subject_dirs{i}, contrast_name);
    end

    desvar(i) = dsn.desvar;
    nscan(i) = dsn.nscan;
    k_cols(i) = dsn.k;
    if ~isempty(dsn.erdf), erdf(i) = dsn.erdf; end
    if ~isempty(dsn.TR), TRs(i) = dsn.TR; end

    % Contrast image: SPM records the file name; find it in this subject's folder
    con_default = sprintf('con_%04d', k);
    if isfield(SPM.xCon(k), 'Vcon') && ~isempty(SPM.xCon(k).Vcon) && isfield(SPM.xCon(k).Vcon, 'fname')
        con_files{i} = find_image_in_dir(subject_dirs{i}, SPM.xCon(k).Vcon.fname, con_default);
    else
        con_files{i} = find_image_in_dir(subject_dirs{i}, '', con_default);
    end

    if isfield(SPM, 'VResMS') && ~isempty(SPM.VResMS) && isfield(SPM.VResMS, 'fname')
        resms_files{i} = find_image_in_dir(subject_dirs{i}, SPM.VResMS.fname, 'ResMS');
    else
        resms_files{i} = find_image_in_dir(subject_dirs{i}, '', 'ResMS');
    end

    if isempty(con_files{i})
        error('canlab_effect_size_map:NoConImage', 'Contrast image for contrast %d not found in %s. Has the contrast been estimated?', k, subject_dirs{i});
    end
end

if verbose
    fprintf('Contrast: %s\n', contrast_name);
    fprintf('Loading %d contrast images\n', N);
end
con_obj = fmri_data(char(con_files), 'noverbose');

wh_no_resms = cellfun(@isempty, resms_files);
if any(wh_no_resms)
    warning('canlab_effect_size_map:NoResMS', 'ResMS images missing for %d subject(s); skipping within/between decomposition.', sum(wh_no_resms));
    resms_obj = [];
else
    if verbose, fprintf('Loading %d ResMS images\n', N); end
    resms_obj = fmri_data(char(resms_files), 'noverbose');
end

design = struct();
design.desvar = desvar;
design.nscan = nscan;
design.k = round(mean(k_cols));
design.erdf = erdf;
if all(isnan(erdf)), design.erdf = []; end
design.TR = mean(TRs, 'omitnan');
if isnan(design.TR), design.TR = []; end
design.contrast_name = contrast_name;
design.contrast_index = k;
design.c = [];
design.source = 'SPM folders';
design.subject_dirs = subject_dirs;

subject_names = subject_dirs;
end


function fname = find_image_in_dir(d, recorded_name, default_base)
% Return the full path of an image in folder d, trying the recorded file
% name first (basename only, in case the analysis folder was moved), then
% default_base with .nii / .img extensions. Returns '' if not found.
fname = '';
candidates = {};
if ~isempty(recorded_name)
    [~, nm, ext] = fileparts(strtok(recorded_name, ','));  % strip ",1" volume suffix
    candidates{end + 1} = [nm ext];
    candidates{end + 1} = [nm '.nii'];
    candidates{end + 1} = [nm '.img'];
end
candidates{end + 1} = [default_base '.nii'];
candidates{end + 1} = [default_base '.img'];

for i = 1:numel(candidates)
    f = fullfile(d, candidates{i});
    if exist(f, 'file') == 2
        fname = f;
        return
    end
end
end



% =========================================================================
% Subfunctions: statistics
% =========================================================================

function pw = power_from_ncp(ncp, df, u, method)
% Power of a one-tailed t-test with critical value u, df degrees of freedom,
% and noncentrality parameter ncp (expected t). ncp may be a vector.
switch method
    case 'noncentral'
        pw = 1 - nctcdf(u, df, ncp);
    case 'shift'
        pw = tcdf(u - ncp, df, 'upper');
    otherwise
        error('canlab_effect_size_map:BadMethod', 'Unknown power_method %s', method);
end
pw(~isfinite(ncp)) = NaN;
end


function u = corrected_threshold(alpha_c, df, nvox, smoothness)
% FWE-corrected critical t (one-tailed): GRF, Bonferroni, or min of both
u_bonf = tinv(1 - alpha_c / nvox, df);

switch smoothness.method
    case 'bonferroni'
        u = u_bonf;
    case 'grf'
        u = spm_uc_RF(alpha_c, [1 df], 'T', smoothness.resels, 1);
    case 'auto'   % smaller of the GRF and Bonferroni thresholds, as SPM does
        u_grf = spm_uc_RF(alpha_c, [1 df], 'T', smoothness.resels, 1);
        u = min(u_grf, u_bonf);
    otherwise
        error('canlab_effect_size_map:BadCorrection', 'Unknown correction method %s', smoothness.method);
end
end


function smoothness = get_smoothness(con, template, fwhm_in, resels_in, correction, verbose)
% Determine resel counts for random field theory thresholds. Estimates
% smoothness from the standardized residuals of the group one-sample
% t-test unless FWHM or resels are supplied.

smoothness = struct('fwhm', [], 'resels', [], 'method', '', 'description', '', 'source', '', 'NISC', [], 'NISC_by_N', []);

have_spm = exist('spm_uc_RF', 'file') == 2 && exist('spm_est_smoothness', 'file') == 2 && exist('spm_resels', 'file') == 2;

if strcmp(correction, 'bonferroni')
    smoothness.method = 'bonferroni';
    smoothness.description = 'Bonferroni over in-mask voxels';
    return
end

if ~have_spm
    if strcmp(correction, 'grf')
        error('canlab_effect_size_map:NoSPM', 'GRF correction requires SPM (spm_uc_RF, spm_est_smoothness, spm_resels) on the path.');
    end
    warning('canlab_effect_size_map:NoSPM', 'SPM not found on the path; using Bonferroni correction.');
    smoothness.method = 'bonferroni';
    smoothness.description = 'Bonferroni over in-mask voxels (SPM not found)';
    return
end

if ~isempty(resels_in)
    R = double(resels_in(:)');
    if isscalar(R), R = [0 0 0 R]; end
    smoothness.resels = R;
    smoothness.fwhm = double(fwhm_in(:)');
    src = 'supplied resel counts';

else
    % Work in a temporary folder: spm_est_smoothness writes RPV.nii to pwd
    tmpdir = tempname;
    mkdir(tmpdir);
    startdir = pwd;
    cleanup = onCleanup(@() cleanup_tmpdir(startdir, tmpdir)); %#ok<NASGU>
    cd(tmpdir);

    % Mask image covering the analyzed voxels
    maskobj = template;
    maskobj.dat = ones(size(template.dat, 1), 1);
    maskfile = fullfile(tmpdir, 'mask.nii');
    write_map(maskobj, maskfile, false);

    if ~isempty(fwhm_in)
        smoothness.fwhm = double(fwhm_in(:)');
        src = 'supplied FWHM';
    else
        % Standardized residuals of the one-sample t-test, following spm_spm
        N = size(con, 2);
        resid = con - mean(con, 2);
        ResSS = sum(resid .^ 2, 2);
        sres = resid ./ sqrt(ResSS ./ (N - 1));
        sres(~isfinite(sres)) = 0;

        nSres = min(N, 64);                      % SPM default: at most 64 residual images
        iRes = round(linspace(1, N, nSres));
        resfiles = cell(nSres, 1);
        for j = 1:nSres
            r = template;
            r.dat = sres(:, iRes(j));
            resfiles{j} = fullfile(tmpdir, sprintf('sres_%03d.nii', j));
            write_map(r, resfiles{j}, false);
        end

        try
            [fwhm, ~, R] = spm_est_smoothness(char(resfiles), maskfile, [N N - 1]);
            smoothness.fwhm = fwhm(:)';
            smoothness.resels = R(:)';
            src = 'estimated from group residuals';
        catch err
            warning('canlab_effect_size_map:SmoothnessFailed', 'Smoothness estimation failed (%s); using Bonferroni correction.', err.message);
            smoothness.method = 'bonferroni';
            smoothness.description = 'Bonferroni over in-mask voxels (smoothness estimation failed)';
            return
        end
    end

    if isempty(smoothness.resels)
        smoothness.resels = spm_resels(smoothness.fwhm, spm_vol(maskfile), 'I');
        smoothness.resels = smoothness.resels(:)';
    end
end

switch correction
    case 'grf'
        smoothness.method = 'grf';
        smoothness.description = 'Random field theory (spm_uc_RF)';
    otherwise
        smoothness.method = 'auto';
        smoothness.description = 'Smaller of random field theory (spm_uc_RF) and Bonferroni thresholds';
end
smoothness.source = src;

if verbose
    if isempty(smoothness.fwhm)
        fwhm_str = 'FWHM not given';
    else
        fwhm_str = sprintf('FWHM = [%3.1f %3.1f %3.1f] voxels', smoothness.fwhm);
    end
    fprintf('Smoothness (%s): %s, %3.1f resels\n', src, fwhm_str, smoothness.resels(end));
end
end


function cleanup_tmpdir(startdir, tmpdir)
cd(startdir);
if exist(tmpdir, 'dir')
    try
        rmdir(tmpdir, 's');
    catch
    end
end
end



% =========================================================================
% Subfunctions: maps, writing, plotting
% =========================================================================

function m = make_map(template, vals, name)
% Create an fmri_data map object from the template and a voxel vector
m = template;
m.dat = single(vals(:));
m.image_names = name;
m.history = {sprintf('canlab_effect_size_map: %s', name)};
end


function write_map(m, fname, verbose)
% Write a map object to disk, quietly unless verbose
m.fullpath = fname;
if verbose
    write(m, 'overwrite');
else
    evalc('write(m, ''overwrite'')');
end
end


function plot_maps(maps, have_within)
% Orthviews of the variance components and power map
try
    if have_within
        obj = maps.residual_std;
        obj.dat = [maps.residual_std.dat maps.within_subjects_std.dat maps.between_subjects_std.dat maps.power_best.dat];
        names = {'Residual std' 'Within-Ss std' 'Between-Ss std' 'FWE corrected power'};
    else
        obj = maps.cohens_d;
        obj.dat = [maps.cohens_d.dat maps.power_best.dat];
        names = {'Cohen''s d' 'FWE corrected power'};
    end
    obj.image_names = char(names);

    orthviews(obj);
    for i = 1:numel(names)
        try, spm_orthviews_name_axis(names{i}, i); catch, end
    end
catch err
    warning('canlab_effect_size_map:PlotFailed', 'Could not display orthviews (%s).', err.message);
end
end


function plot_allocation(results, have_within)
% 2 x 2 figure: allocation, variance pie, power curves, power distribution
alloc = results.allocation;
best = alloc.best;

create_figure('Optimal allocation', 2, 2);

% Sessions and functional hours vs N
subplot(2, 2, 1);
plot(alloc.N, alloc.sessions_per_subject, 'k', 'LineWidth', 3);
plot(alloc.N, alloc.functional_hours_per_subject, 'k:', 'LineWidth', 3);
legend({'Sessions per subject' 'Functional scan hours per subject'}, 'Location', 'best');
axis tight
plot_vertical_line(alloc.session_change_N);
xlabel('Number of subjects');
title(sprintf('Allocation of %3.0f scanner hours', alloc.hours));

% Pie chart of within vs. between variance at the best allocation
subplot(2, 2, 2);
if have_within
    hh = pie([mean(best.sig2_within_at_best, 'omitnan') mean(results.sig2_between, 'omitnan')], {'Within', 'Between'});
    hp = findobj(hh, 'Type', 'Patch');
    if numel(hp) >= 2
        set(hp(1), 'FaceColor', [.3 .3 .3]);
        set(hp(2), 'FaceColor', [.8 .8 .8]);
    end
    title(sprintf('Variance at N = %d', best.N));
else
    text(.1, .5, 'No within-subject variance information', 'FontSize', 12);
    axis off
end

% Power curves
subplot(2, 2, 3);
plot(alloc.N, alloc.mean_power_corrected, 'k', 'LineWidth', 3);
legstr = {'FWE corrected'};
grays = linspace(.45, .75, numel(alloc.alpha_uncorrected));
for a = 1:numel(alloc.alpha_uncorrected)
    plot(alloc.N, alloc.mean_power_uncorrected(:, a), '-', 'Color', grays(a) * [1 1 1], 'LineWidth', 3);
    legstr{end + 1} = sprintf('p < %g', alloc.alpha_uncorrected(a)); %#ok<AGROW>
end
xlabel('Number of subjects');
ylabel('Mean power in search area');
legend(legstr, 'Location', 'best');
axis tight
plot_vertical_line(alloc.session_change_N);
plot_vertical_line(best.N, 'r');
hh = plot_horizontal_line(.8); set(hh, 'LineStyle', ':');
title(sprintf('Best: N = %d, %d min/subject', best.N, best.functional_min_per_subject));

% Distribution of corrected power across voxels at the best allocation
subplot(2, 2, 4);
histogram(results.power_map, 20, 'FaceColor', [.5 .5 .5]);
xlabel('FWE-corrected power');
ylabel('Voxels');
title('Power across search area');
end
