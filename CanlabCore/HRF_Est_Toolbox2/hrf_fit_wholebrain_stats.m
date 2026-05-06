function out = hrf_fit_wholebrain_stats(fmri_nii, events_tsv, varargin)
%HRF_FIT_WHOLEBRAIN_STATS Vectorized whole-brain FIR/sFIR HRF beta and T maps.
%
% out = hrf_fit_wholebrain_stats(fmri_nii, events_tsv, ...)
%
% This fits an FIR-style time-unfolding model to every voxel in a 4D fMRI
% image and returns two CANlab statistic_image objects:
%   out.b  - beta/HRF amplitude maps, one 3D volume per condition x lag
%   out.t  - T maps for the same condition x lag volumes
%
% Both objects include .ste, .p, .dfe, .N, .image_labels, and .volInfo.
% If OutputPrefix is provided, 4D NIfTI files are written to disk.

p = inputParser;
p.addRequired('fmri_nii', @(x) ischar(x) || isstring(x) || isa(x, 'fmri_data'));
p.addRequired('events_tsv', @(x) ischar(x) || isstring(x) || istable(x));
p.addParameter('TR', [], @(x) isempty(x) || (isscalar(x) && x > 0));
p.addParameter('MaskNii', '', @(x) ischar(x) || isstring(x) || isa(x, 'fmri_mask_image'));
p.addParameter('Conditions', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('WindowSeconds', 30, @(x) isnumeric(x) && all(x(:) > 0));
p.addParameter('Mode', 'FIR', @(x) ischar(x) || isstring(x));
p.addParameter('Nuisance', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('OutputPrefix', '', @(x) ischar(x) || isstring(x));
p.addParameter('Overwrite', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('PThresh', [], @(x) isempty(x) || (isscalar(x) && x > 0 && x < 1));
p.addParameter('ThreshType', 'unc', @(x) ischar(x) || isstring(x));
p.addParameter('WriteThresholdedT', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('ChunkSize', 50000, @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x, 1) == 0);
p.addParameter('ScaleMode', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(fmri_nii, events_tsv, varargin{:});
opts = p.Results;

if isa(fmri_nii, 'fmri_data')
    data_obj = fmri_nii;
    fmri_name = '';
    if isempty(opts.TR)
        if isprop(data_obj, 'image_metadata') && isfield(data_obj.image_metadata, 'TR_in_sec') && ...
                ~isnan(data_obj.image_metadata.TR_in_sec)
            TR = data_obj.image_metadata.TR_in_sec;
        else
            error('TR is required when fmri_nii is an fmri_data object without image_metadata.TR_in_sec.');
        end
    else
        TR = opts.TR;
    end
    n_tp = size(data_obj.dat, 2);
else
    fmri_name = char(fmri_nii);
    if ~exist(fmri_name, 'file'), error('fMRI file not found: %s', fmri_name); end
    [TR, n_tp] = local_get_tr_and_ntp(fmri_name, opts.TR);

    mask_input = opts.MaskNii;
    if isstring(mask_input), mask_input = char(mask_input); end

    if isempty(mask_input)
        data_obj = fmri_data(fmri_name, 'noverbose');
    else
        data_obj = fmri_data(fmri_name, mask_input, 'noverbose');
    end
end

if istable(events_tsv)
    E = events_tsv;
else
    events_name = char(events_tsv);
    if ~exist(events_name, 'file'), error('Events file not found: %s', events_name); end
    E = hrf_load_events_tsv(events_name);
end

if isempty(opts.Conditions)
    cond_names = unique(E.trial_type, 'stable');
else
    cond_names = cellstr(string(opts.Conditions));
end

Runc = hrf_build_stick_functions(E, cond_names, TR, n_tp);
[X, design_info] = local_build_fir_design(Runc, cond_names, TR, opts.WindowSeconds, opts.Mode, opts.Nuisance);

if size(X, 1) ~= size(data_obj.dat, 2)
    error('Design rows (%d) must match fMRI time points (%d).', size(X, 1), size(data_obj.dat, 2));
end

[PX, pen] = local_get_pseudoinverse(X, design_info, opts.Mode);
hat_trace = trace(X * PX);
dfe = max(size(X, 1) - hat_trace, 1);
coef_var_scale = max(diag(PX * PX'), 0);

keep = design_info.fir_columns;
n_keep = numel(keep);
n_vox = size(data_obj.dat, 1);

beta_dat = zeros(n_vox, n_keep, 'single');
ste_dat = zeros(n_vox, n_keep, 'single');
t_dat = zeros(n_vox, n_keep, 'single');
p_dat = zeros(n_vox, n_keep, 'single');

chunk_size = min(opts.ChunkSize, n_vox);
if opts.Verbose
    fprintf('Fitting whole-brain %s model: %d voxels, %d time points, %d output maps\n', ...
        upper(char(opts.Mode)), n_vox, size(X, 1), n_keep);
end

for first_vox = 1:chunk_size:n_vox
    last_vox = min(first_vox + chunk_size - 1, n_vox);
    wh = first_vox:last_vox;
    Y = double(data_obj.dat(wh, :)');

    switch lower(char(opts.ScaleMode))
        case 'none'
            % leave data in native units
        case {'zscore', 'z'}
            Y = zscore(Y, 0, 1);
            Y(isnan(Y)) = 0;
        otherwise
            error('Unknown ScaleMode: %s. Use ''none'' or ''zscore''.', char(opts.ScaleMode));
    end

    B = PX * Y;
    R = Y - X * B;
    mse = sum(R .^ 2, 1) ./ dfe;
    SE = sqrt(coef_var_scale * mse);

    Bkeep = B(keep, :);
    SEkeep = SE(keep, :);
    Tkeep = Bkeep ./ SEkeep;
    Pkeep = 2 * (1 - tcdf(abs(Tkeep), dfe));
    Pkeep(Pkeep == 0) = eps;

    beta_dat(wh, :) = single(Bkeep');
    ste_dat(wh, :) = single(SEkeep');
    t_dat(wh, :) = single(Tkeep');
    p_dat(wh, :) = single(Pkeep');

    if opts.Verbose
        fprintf('  voxels %d-%d / %d\n', first_vox, last_vox, n_vox);
    end
end

labels = design_info.labels;
meta_table = design_info.metadata_table;

b_obj = statistic_image;
b_obj.type = sprintf('%s HRF beta', upper(char(opts.Mode)));
b_obj.dat = beta_dat;
b_obj.p = p_dat;
b_obj.p_type = sprintf('Two-tailed P-values from beta/SE, dfe = %.3f', dfe);
b_obj.ste = ste_dat;
b_obj.sig = true(size(beta_dat));
b_obj.N = size(X, 1);
b_obj.dfe = dfe;
b_obj.volInfo = data_obj.volInfo;
b_obj.removed_voxels = data_obj.removed_voxels;
b_obj.removed_images = false(1, n_keep);
b_obj.image_labels = labels;
b_obj.dat_descrip = sprintf('Whole-brain %s HRF beta values: condition x lag maps', upper(char(opts.Mode)));

t_obj = statistic_image;
t_obj.type = 'T';
t_obj.dat = t_dat;
t_obj.p = p_dat;
t_obj.p_type = sprintf('Two-tailed P-values from beta/SE, dfe = %.3f', dfe);
t_obj.ste = ste_dat;
t_obj.sig = true(size(t_dat));
t_obj.N = size(X, 1);
t_obj.dfe = dfe;
t_obj.volInfo = data_obj.volInfo;
t_obj.removed_voxels = data_obj.removed_voxels;
t_obj.removed_images = false(1, n_keep);
t_obj.image_labels = labels;
t_obj.dat_descrip = sprintf('Whole-brain %s HRF T values: condition x lag maps', upper(char(opts.Mode)));

if ~isempty(opts.PThresh)
    t_obj = threshold(t_obj, opts.PThresh, char(opts.ThreshType), 'noverbose');
end

out = struct();
out.b = b_obj;
out.t = t_obj;
out.design_matrix = X;
out.design_info = design_info;
out.design_info.penalty = pen;
out.metadata_table = meta_table;
out.conditions = cond_names;
out.TR = TR;
out.dfe = dfe;
out.N = size(X, 1);
out.input_fmri = fmri_name;
out.paths = struct();

if ~isempty(opts.OutputPrefix)
    out.paths = local_write_outputs(out, char(opts.OutputPrefix), logical(opts.Overwrite), logical(opts.WriteThresholdedT));
end
end

function [TR, n_tp] = local_get_tr_and_ntp(fmri_name, requested_tr)
info = niftiinfo(fmri_name);
sz = info.ImageSize;
if numel(sz) < 4
    error('Expected a 4D fMRI image, got %d dimensions.', numel(sz));
end
n_tp = sz(4);

TR = requested_tr;
if isempty(TR) && isfield(info, 'PixelDimensions') && numel(info.PixelDimensions) >= 4
    TR = info.PixelDimensions(4);
end
if isempty(TR) || TR <= 0
    error('Could not infer TR from NIfTI header. Pass ''TR'' explicitly.');
end
end

function [X, info] = local_build_fir_design(Runc, cond_names, TR, window_seconds, mode, nuisance)
numstim = numel(Runc);
len = numel(Runc{1});
if isscalar(window_seconds)
    window_seconds = repmat(window_seconds, 1, numstim);
end
if numel(window_seconds) ~= numstim
    error('WindowSeconds must be scalar or one value per condition.');
end

Runs = zeros(len, numstim);
for i = 1:numstim
    Runs(:, i) = Runc{i}(:);
end

DX_all = cell(1, numstim);
tlen_all = zeros(1, numstim);
fir_columns_by_condition = cell(1, numstim);
labels = {};
condition_col = {};
condition_index = [];
lag_index = [];
lag_seconds = [];

start_col = 1;
for i = 1:numstim
    t = 1:TR:window_seconds(i);
    tlen_all(i) = numel(t);
    DX_i = tor_make_deconv_mtx3(Runs(:, i), tlen_all(i), 1);
    DX_all{i} = DX_i(:, 1:tlen_all(i));

    fir_columns_by_condition{i} = start_col:(start_col + tlen_all(i) - 1);
    for j = 1:tlen_all(i)
        label = sprintf('%s_lag%03d_%0.3gs', local_safe_label(cond_names{i}), j, t(j));
        labels{end + 1, 1} = label; %#ok<AGROW>
        condition_col{end + 1, 1} = cond_names{i}; %#ok<AGROW>
        condition_index(end + 1, 1) = i; %#ok<AGROW>
        lag_index(end + 1, 1) = j; %#ok<AGROW>
        lag_seconds(end + 1, 1) = t(j); %#ok<AGROW>
    end
    start_col = start_col + tlen_all(i);
end

DX = horzcat(DX_all{:});
X = [DX ones(len, 1)];

if ~isempty(nuisance)
    if size(nuisance, 1) ~= len
        error('Nuisance matrix must have %d rows to match fMRI time points.', len);
    end
    X = [X nuisance];
end

info = struct();
info.mode = char(mode);
info.TR = TR;
info.tlen = tlen_all;
info.fir_columns = 1:sum(tlen_all);
info.fir_columns_by_condition = fir_columns_by_condition;
info.intercept_column = sum(tlen_all) + 1;
info.nuisance_columns = (info.intercept_column + 1):size(X, 2);
info.labels = labels;
info.metadata_table = table((1:numel(labels))', condition_col, condition_index, lag_index, lag_seconds, labels, ...
    'VariableNames', {'volume_index', 'condition', 'condition_index', 'lag_index', 'lag_seconds', 'image_label'});
end

function [PX, pen] = local_get_pseudoinverse(X, info, mode)
mode = lower(char(mode));
pen = zeros(size(X, 2));

switch mode
    case 'fir'
        PX = pinv(X);
    case 'sfir'
        start_idx = 1;
        for i = 1:numel(info.tlen)
            tlen = info.tlen(i);
            C = (1:tlen)' * ones(1, tlen);
            h = sqrt(1 / (7 / info.TR));
            v = 0.1;
            sig = 1;
            R = v * exp(-h / 2 * (C - C') .^ 2);
            RI = R \ eye(size(R));
            end_idx = start_idx + tlen - 1;
            pen(start_idx:end_idx, start_idx:end_idx) = sig ^ 2 * RI;
            start_idx = end_idx + 1;
        end
        PX = (X' * X + pen) \ X';
    otherwise
        error('Unknown Mode: %s. Use ''FIR'' or ''sFIR''.', mode);
end
end

function paths = local_write_outputs(out, prefix, overwrite, write_thresholded_t)
paths = struct();
paths.beta = [prefix '_beta.nii'];
paths.t = [prefix '_t.nii'];
paths.metadata = [prefix '_metadata.csv'];

write_args = {'fname', paths.beta};
if overwrite, write_args{end + 1} = 'overwrite'; end
write(out.b, write_args{:});

write_args = {'fname', paths.t};
if overwrite, write_args{end + 1} = 'overwrite'; end
write(out.t, write_args{:});

writetable(out.metadata_table, paths.metadata);

if write_thresholded_t
    paths.t_thresholded = [prefix '_t_thresh.nii'];
    write_args = {'fname', paths.t_thresholded, 'thresh'};
    if overwrite, write_args{end + 1} = 'overwrite'; end
    write(out.t, write_args{:});
end
end

function s = local_safe_label(s)
s = char(s);
s = matlab.lang.makeValidName(s);
end
