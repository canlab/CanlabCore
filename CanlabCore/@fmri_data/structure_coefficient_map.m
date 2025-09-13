function [r1_obj, r2_obj, rdiff_obj] = structure_coefficient_map(obj, pattern1, varargin)
% Calculate structure coefficient map (aka model encoding map) for linear model pattern1 using data in obj
%
% [r1_obj, r2_obj, rdiff_obj] = structure_coefficient_map(obj, pattern1, varargin)
%
% pattern1 is an fmri_data object with model weights for a linear multivariate model
%
% obj is an fmri_data object with a dataset to calculate structure
% coefficients. For valid statistics, images should be independent (e.g., 1
% per participant).
%
% Optional inputs:
% 'pattern2' followed by 2nd model to compare
% 'montage'     Plot montage(s) of results
% 'doplots'     Plot montage(s) of results
% 'nboot'       followed by num bootstrap samples, default 100 (just for
% testing purposes). Recommended: 5000 or more for final analysis
%
% Example: Single model
% ---------------------------------------------------------
% ncs = load_image_set('ncs');
% obj = load(which('anic_2021_all_current_v20220831.mat'));
% obj = load_image_set('emotionreg');
% r1_obj = structure_coefficient_map(obj, ncs);
%
%
% Example: Comparison of two models
% ---------------------------------------------------------
% nps = load_image_set('nps');
% vps = load_image_set('vps');
% obj = load_image_set('emotionreg');
%
% pattern1 = nps;  pattern2 = vps;
% [ncs_obj, nps_obj, diff_obj] = structure_coefficient_map(obj, ncs, 'pattern2', nps, 'montage', 'nboot', 1000);
%
% % Re-threshold maps and display again:
% r1_obj = threshold(r1_obj, .05, 'unc');
% rdiff_obj = threshold(rdiff_obj, .05, 'unc');
%

% SEE ALSO
% anic_2021_all_current_v20220831.mat
% get_model_encoding_map.m
% anic_model_encoding_mpa2.m
% prep_3_create_model_encoding_maps.m
% create_model_encoding_maps.m
% method:
% bootstrap_structure_coeff_diff.m
%
% notes structure coeffs source reconstruction.rtfd

% ----------------------------------------------------------------------
% Parse inputs
% ----------------------------------------------------------------------
% This 2019 version uses the inputParser object. Older schemes are below.
% Note: With this, you can pass in EITHER keyword, value pairs OR a
% structure with fields that define keywords.
% Note: You can also make this a subfunction at the end of your function,
% e.g., function ARGS = parse_inputs(varargin)...end, and call it using:

% Keywords - these should be removed before parsing
% ----------------------------------------------------------------------
wh = strcmp(varargin, 'montage');
if any(strcmp(varargin, 'montage'))
    varargin(wh) = [];
    varargin{end + 1} = 'doplots';
    varargin{end + 1} = 1;
end

ARGS = parse_inputs(varargin{:});

% If you want to distribute arguments back out to variables, use this:

fn = fieldnames(ARGS);

for i = 1:length(fn)
    str = sprintf('%s = ARGS.(''%s'');', fn{i}, fn{i});
    eval(str)
end

% Logical flags - these override earlier inputs
% This doesn't work tho
% ----------------------------------------------------------------------
% if any(strcmp(varargin, 'noverbose')), doverbose = false; end
% if any(strcmp(varargin, 'noplots')), doplots = false; end
% if any(strcmp(varargin, 'noboot')), doboot = false; end

% ----------------------------------------------------------------------
% Prep
% ----------------------------------------------------------------------

% outputs
% r1_obj = [];
r2_obj = [];
rdiff_obj = [];

obj.dat = double(obj.dat);  % superstition

if size(pattern1.dat, 2) > 1
    error('Pattern 1 (model) must have only 1 image')
end

if ~isempty(pattern2) && size(pattern1.dat, 2) > 1
    error('Pattern 2 (model) must have only 1 image')
end


% ----------------------------------------------------------------------
% Apply patterns
% ----------------------------------------------------------------------

pexp1 = double(apply_mask(obj, pattern1, 'pattern_expression', 'cosine_similarity'));

% Normalize - z-score so we are computing correlations (or partial
% correlations)
pexp1 = zscore(pexp1);

if ~isempty(pattern2)

    pexp2 = double(apply_mask(obj, pattern2, 'pattern_expression', 'cosine_similarity'));

    pexp2 = zscore(pexp2);

end


% ----------------------------------------------------------------------
% Calculate structure coefficient maps for full sample
% ----------------------------------------------------------------------

r1 = corr(obj.dat', pexp1);

if ~isempty(pattern2)

    r2 = corr(obj.dat', pexp2);

    rdiff = r1 - r2;

end


% ----------------------------------------------------------------------
% Bootstrapping
% ----------------------------------------------------------------------

%% Set up bootstrap samples

bootsam = setup_boot_samples(pexp1, nboot);

%%

nvox = size(obj.dat, 1);
[r1_boot, r2_boot, rdiff_boot] = deal(NaN .* zeros(nvox, nboot));

tic
fprintf('Bootstrapping, %3.0f samples\n', nboot)

for i = 1:nboot

    indx = bootsam(:, i);

    obj_boot = get_wh_image(obj, indx);

    pexp1_boot = pexp1(indx);

    if ~isempty(pattern2)

        pexp2_boot = pexp2(indx);

    end

    % Get summary statistic maps for this bootstrap sample
    % voxels x bootstrap samples, each column is a map
    r1_boot(:, i) = corr(double(obj_boot.dat'), pexp1_boot);

    if ~isempty(pattern2)

        r2_boot(:, i) = corr(double(obj_boot.dat'), pexp2_boot);

        rdiff_boot(:, i) = r1_boot(:, i) - r2_boot(:, i);

    end

end % bootstrap iteration

toc

% ----------------------------------------------------------------------
% Reconstruct statistic images
% ----------------------------------------------------------------------


r1_obj = get_statistic_image(r1_boot, r1, obj, 'Structure coeff map for pattern 1');

if ~isempty(pattern2)

    r2_obj = get_statistic_image(r2_boot, r2, obj, 'Structure coeff map for pattern 2');

    rdiff_obj = get_statistic_image(rdiff_boot, rdiff, obj, 'Structure coeff map for pattern 1 - 2');

end


% ----------------------------------------------------------------------
% Plots
% ----------------------------------------------------------------------


if doplots

    % single
%     create_figure('montage1'); axis off;
    mylabels = {'Pattern 1, unthresholded' 'Thresholded, q < 0.05 FDR'};
    show_montage(r1_obj, mylabels);

    if ~isempty(pattern2)

%         create_figure('montage2'); axis off;
        mylabels = {'Pattern 2, unthresholded' 'Thresholded, q < 0.05 FDR'};
        show_montage(r2_obj, mylabels);

%         create_figure('montagediff'); axis off;
        mylabels = {'Pattern 1-2, unthresholded' 'Thresholded, q < 0.05 FDR'};
        show_montage(rdiff_obj, mylabels);

    end


end


end % main function



%% Subfunctions

function r1_obj = get_statistic_image(r1_boot, r1, obj, descrip)

r1_p = bootbca_pval(0, @mean, r1_boot' ,r1', obj.dat');
r1_p = r1_p';
% r1_z = r1_z';

r1_fdr_thresh = FDR(r1_p, 0.05);
if isempty(r1_fdr_thresh),  r1_fdr_thresh = -Inf; end

r1_obj = statistic_image('dat', r1, ...
    'volInfo', obj.volInfo, ...
    'p', r1_p, ...
    'sig', r1_p < r1_fdr_thresh, ...
    'ste', [], ...
    'dat_descrip', descrip, ...
    'removed_voxels', obj.removed_voxels);

end




function show_montage(r1_obj, mylabels)
r1_to_plot = r1_obj;
r1_to_plot.p(:, 2) = r1_obj.p;
r1_to_plot.sig(:, 2) = r1_obj.sig;
r1_to_plot.dat(:, 2) = r1_obj.dat;
r1_to_plot.sig(:, 1) = 1; % unthresholded
r1_to_plot.image_labels = mylabels;
montage(r1_to_plot);
end





function ARGS = parse_inputs(varargin)

p = inputParser;

% Validation functions - customized for each type of input
% ----------------------------------------------------------------------

valfcn_scalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar'});

% valfcn_number = @(x) validateattributes(x, {'numeric'}, {'nonempty'}); % scalar or vector

% valfcn_fmridata = @(x) isempty(x) || isa(x, 'fmri_data') ;
valfcn_imagevec = @(x) isempty(x) || isa(x, 'image_vector') ;

% Validation: [x1 x2 x3] triplet
valfcn_xyz = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'size', [1 3]});

valfcn_logical = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', '>=', 0, '<=', 1}); % could enter numeric 0,1 or logical


% Required inputs
% ----------------------------------------------------------------------
% p.addRequired('x', valfcn_custom);

% Optional inputs
% ----------------------------------------------------------------------
% Pattern: keyword, value, validation function handle

p.addParameter('color', [.9 .2 0], valfcn_xyz);
p.addParameter('pattern2', [], valfcn_imagevec); % can be scalar or vector
p.addParameter('doplots', 0, valfcn_logical);
p.addParameter('doboot', 0, valfcn_logical);
p.addParameter('nboot', 100, valfcn_scalar);

% Parse inputs and distribute out to variable names in workspace
% ----------------------------------------------------------------------
p.parse(varargin{:});

ARGS = p.Results;

end % parse_inputs

