function [obj_out, statstab] = normalize_gm_shift_scale(obj, varargin)
% Normalize gray-matter voxel intensities across subjects in an fmri_data
% object by:
%   (1) Removing a subject-specific additive shift estimated from CSF/WM
%       medians
%   (2) Correcting a subject-specific multiplicative scale estimated from
%       MADs in GM, WM, and CSF
%
% The method:
%   - Uses canonical masks for GM, WM, and CSF:
%       {'gray_matter_mask_sparse.img', ...
%        'canonical_white_matter.img', ...
%        'canonical_ventricles.img'}
%   - Resamples masks into the space of the fmri_data object using
%     RESAMPLE_SPACE to ensure voxel alignment.
%   - Calls NORMALIZE_GM_SHIFT_SCALE (voxel-level function) on the data.
%   - Returns:
%       * an fmri_data object with GM voxels shift- and scale-normalized,
%         and non-GM voxels left unchanged except that they are not further
%         processed
%       * a MATLAB table containing all values from the STATS output of
%         NORMALIZE_GM_SHIFT_SCALE, appended to the existing
%         metadata_table.
%
% USAGE
%   [obj_out, statstab] = gm_shift_scale_normalize(obj, ...
%                                'log_scale', false, 'trim_pct', 5, ...
%                                'mask_files', masks_cell);
%
% INPUTS
%   obj       : fmri_data object with .dat of size [V x S]
%               - V = number of voxels
%               - S = number of images/subjects
%
% OPTIONAL NAME/VALUE INPUTS
%   'log_scale' : logical (default = false)
%                 - Passed to NORMALIZE_GM_SHIFT_SCALE.
%                 - If true: log-scale regression for scale model
%                   log(r_GM) ~ log(r_CSF) + log(r_WM)
%                 - If false: linear scale regression
%                   r_GM ~ r_CSF + r_WM
%
%   'trim_pct'  : scalar (default = 5)
%                 - Percentage trimmed from lower and upper tails when
%                   estimating medians/MADs in each tissue.
%
%   'mask_files': 1 x 3 cell array of mask filenames
%                 (default):
%                 {'gray_matter_mask_sparse.img', ...
%                  'canonical_white_matter.img', ...
%                  'canonical_ventricles.img'}
%
% OUTPUTS
%   obj_out  : fmri_data object
%              - .dat is V x S, with GM voxels normalized by the shift/scale
%                model; non-GM voxels are copied from the input.
%              - .metadata_table is the input metadata_table with new
%                columns appended containing all (subject-level) values from
%                STATS.
%
%   statstab : MATLAB table (S rows)
%              - All values from STATS struct are represented as columns.
%              - This same table is appended to obj_out.metadata_table.
%
% NOTES
%   - This method does not reduce the number of voxels; it keeps the full
%     spatial geometry but only modifies GM voxels in .dat.
%   - If you prefer to hard-mask to GM (retain only GM voxels in .dat),
%     you can add an additional step to subset voxels by the GM mask.
%
%   Author:  (adapted for CANlab-style documentation)
%   Date:    2025-12-09
%

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

p = inputParser;
p.addRequired('obj', @(x) isa(x, 'fmri_data'));

default_masks = {'gray_matter_mask_sparse.img', ...
                 'canonical_white_matter.img', ...
                 'canonical_ventricles.img'};

p.addParameter('log_scale', false, @(x) islogical(x) && isscalar(x));
p.addParameter('trim_pct',  5,     @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 50);
p.addParameter('mask_files', default_masks, @(x) iscell(x) && numel(x) == 3);

p.parse(obj, varargin{:});

log_scale  = p.Results.log_scale;
trim_pct   = p.Results.trim_pct;
mask_files = p.Results.mask_files;

% -------------------------------------------------------------------------
% Basic checks
% -------------------------------------------------------------------------

if isempty(obj.dat)
    error('gm_shift_scale_normalize:EmptyData', ...
        'The fmri_data object has an empty .dat field.');
end

[V, S] = size(obj.dat);

% -------------------------------------------------------------------------
% Load and resample masks to object space
% -------------------------------------------------------------------------

% GM
gm_mask_obj = fmri_data(mask_files{1});
gm_mask_obj = resample_space(gm_mask_obj, obj, 'nearest');
gm_mask = gm_mask_obj.dat ~= 0 & ~isnan(gm_mask_obj.dat);

% WM
wm_mask_obj = fmri_data(mask_files{2});
wm_mask_obj = resample_space(wm_mask_obj, obj, 'nearest');
wm_mask = wm_mask_obj.dat ~= 0 & ~isnan(wm_mask_obj.dat);

% CSF
csf_mask_obj = fmri_data(mask_files{3});
csf_mask_obj = resample_space(csf_mask_obj, obj, 'nearest');
csf_mask = csf_mask_obj.dat ~= 0 & ~isnan(csf_mask_obj.dat);

% Make column vectors and check length
gm_mask  = gm_mask(:);
wm_mask  = wm_mask(:);
csf_mask = csf_mask(:);

if length(gm_mask) ~= V || length(wm_mask) ~= V || length(csf_mask) ~= V
    error('gm_shift_scale_normalize:MaskSizeMismatch', ...
        'Resampled masks do not match the number of voxels in obj.dat.');
end

% -------------------------------------------------------------------------
% Call voxel-level normalization function
% -------------------------------------------------------------------------

Y = obj.dat;  % V x S

[Z, stats] = normalize_gm_shift_scale(Y, gm_mask, wm_mask, csf_mask, ...
                                     'log_scale', log_scale, ...
                                     'trim_pct',  trim_pct);

% -------------------------------------------------------------------------
% Build stats table from STATS struct
% -------------------------------------------------------------------------

statstab = stats_struct_to_table(stats, S);

% -------------------------------------------------------------------------
% Construct output fmri_data object and metadata table
% -------------------------------------------------------------------------

obj_out = obj;            % copy all properties
obj_out.dat = Z;          % normalized data

% Metadata table: append stats columns
if istable(obj.metadata_table)
    mt_in = obj.metadata_table;
    if height(mt_in) ~= S
        warning('gm_shift_scale_normalize:MetadataHeightMismatch', ...
            ['metadata_table height (%d) does not match number of images (%d). ', ...
             'Stats table will still be appended; please verify correspondence.'], ...
             height(mt_in), S);
    end
    obj_out.metadata_table = [mt_in statstab];
else
    % If no metadata_table, create a new one from stats
    obj_out.metadata_table = statstab;
end

end % function gm_shift_scale_normalize

% =========================================================================
% Helper: convert STATS struct to table with S rows
% =========================================================================
function T = stats_struct_to_table(stats, S)
% STATS_STRUCT_TO_TABLE
%
% Convert a STATS struct (from normalize_gm_shift_scale) into a table with
% S rows (one per subject) and appropriately named columns.

fn = fieldnames(stats);
vars = struct();  % accumulate columns here

for k = 1:numel(fn)
    fname = fn{k};
    val   = stats.(fname);

    if isempty(val)
        % Skip empty fields
        continue;
    end

    sz = size(val);

    % Case 1: scalar -> replicate across subjects
    if isscalar(val)
        col = repmat(val, S, 1);
        vars.(fname) = col;

    % Case 2: column vector of length S
    elseif isvector(val) && sz(1) == S && sz(2) == 1
        vars.(fname) = val;

    % Case 3: row vector of length S (unlikely but handle)
    elseif isvector(val) && sz(2) == S && sz(1) == 1
        vars.(fname) = val(:);

    % Case 4: S x K matrix -> create one column per K
    elseif ndims(val) == 2 && sz(1) == S && sz(2) > 1
        K = sz(2);
        for j = 1:K
            newname = sprintf('%s_%d', fname, j);
            vars.(newname) = val(:, j);
        end

    % Case 5: 1 x K vector (e.g., beta_mean, lambda) -> replicate as K columns
    elseif ndims(val) == 2 && sz(1) == 1 && sz(2) > 1
        K = sz(2);
        for j = 1:K
            newname = sprintf('%s_%d', fname, j);
            vars.(newname) = repmat(val(1, j), S, 1);
        end

    else
        % Higher-dimensional or mismatched; skip with a warning
        warning('stats_struct_to_table:SkipField', ...
            'Skipping field "%s" with size [%s].', fname, num2str(sz));
    end

end

T = struct2table(vars);

end