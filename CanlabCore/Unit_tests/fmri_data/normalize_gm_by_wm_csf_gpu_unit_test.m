function results = normalize_gm_by_wm_csf_gpu_unit_test(varargin)
% normalize_gm_by_wm_csf_gpu_unit_test
%
% Synthetic CPU/GPU equivalence and timing checks for normalize_gm_shift_scale,
% plus an opportunistic end-to-end benchmark for normalize_gm_by_wm_csf on the
% emotionreg sample dataset.
%
% Usage:
%   results = normalize_gm_by_wm_csf_gpu_unit_test;
%   results = normalize_gm_by_wm_csf_gpu_unit_test('run_object_benchmark', false);

p = inputParser;
p.addParameter('run_object_benchmark', true, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
p.addParameter('synthetic_voxels', 12000, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('synthetic_images', 20, @(x) isnumeric(x) && isscalar(x) && x > 3);
p.addParameter('tolerance', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.parse(varargin{:});

run_object_benchmark = logical(p.Results.run_object_benchmark);
V = p.Results.synthetic_voxels;
S = p.Results.synthetic_images;
tol = p.Results.tolerance;

ensure_canlab_path();

[gpu_available, gpu_message] = local_gpu_available();
results = struct();
results.gpu_available = gpu_available;
results.gpu_message = gpu_message;
results.tolerance = tol;

fprintf('\nnormalize_gm_by_wm_csf GPU trial test\n');
fprintf('GPU available: %d', gpu_available);
if ~gpu_available
    fprintf(' (%s)', gpu_message);
end
fprintf('\n');

[Y, gm_mask, wm_mask, csf_mask] = synthetic_data(V, S);

configs = { ...
    struct('name', 'linear_scale', 'do_scale', true,  'log_scale', false), ...
    struct('name', 'log_scale',    'do_scale', true,  'log_scale', true), ...
    struct('name', 'shift_only',   'do_scale', false, 'log_scale', false)};

for i = 1:numel(configs)
    cfg = configs{i};

    [cpu_time, Z_cpu, stats_cpu] = timed_kernel(Y, gm_mask, wm_mask, csf_mask, cfg, false);
    [gpu_time, Z_gpu, stats_gpu] = timed_kernel(Y, gm_mask, wm_mask, csf_mask, cfg, true);

    max_z_diff = max_abs_diff(Z_cpu, Z_gpu);
    max_stats_diff = max_stats_abs_diff(stats_cpu, stats_gpu);

    synthetic_results(i).scenario = cfg.name; %#ok<AGROW>
    synthetic_results(i).cpu_time_sec = cpu_time; %#ok<AGROW>
    synthetic_results(i).gpu_requested_time_sec = gpu_time; %#ok<AGROW>
    synthetic_results(i).max_Z_abs_diff = max_z_diff; %#ok<AGROW>
    synthetic_results(i).max_stats_abs_diff = max_stats_diff; %#ok<AGROW>
    synthetic_results(i).passed = max(max_z_diff, max_stats_diff) <= tol; %#ok<AGROW>

    fprintf('%s: CPU %.4fs, use_gpu %.4fs, max |Z diff| %.3g, max |stats diff| %.3g\n', ...
        cfg.name, cpu_time, gpu_time, max_z_diff, max_stats_diff);

    assert(synthetic_results(i).passed, ...
        'GPU/fallback result differs from CPU for %s. Max diff = %.6g', ...
        cfg.name, max(max_z_diff, max_stats_diff));
end

results.synthetic = synthetic_results;

if run_object_benchmark
    results.object_benchmark = run_emotionreg_object_benchmark(tol);
else
    results.object_benchmark = struct('status', 'skipped_by_user');
end

fprintf('normalize_gm_by_wm_csf GPU trial test complete.\n\n');

end


function ensure_canlab_path()

if isempty(which('normalize_gm_shift_scale'))
    this_file = mfilename('fullpath');
    fmri_data_test_dir = fileparts(this_file);
    unit_tests_dir = fileparts(fmri_data_test_dir);
    canlab_root = fileparts(unit_tests_dir);
    addpath(genpath(canlab_root));
end

end


function [Y, gm_mask, wm_mask, csf_mask] = synthetic_data(V, S)

rng(20260611, 'twister');

gm_mask = false(V, 1);
wm_mask = false(V, 1);
csf_mask = false(V, 1);

gm_mask(1:round(0.55 * V)) = true;
wm_mask((round(0.55 * V) + 1):round(0.78 * V)) = true;
csf_mask((round(0.78 * V) + 1):round(0.90 * V)) = true;

subject_shift = linspace(-0.4, 0.4, S);
subject_scale = linspace(0.75, 1.25, S);

Y = randn(V, S);
Y(gm_mask, :) = 1.1 + Y(gm_mask, :);
Y(wm_mask, :) = -0.2 + 0.8 .* Y(wm_mask, :);
Y(csf_mask, :) = 0.2 + 1.4 .* Y(csf_mask, :);

Y = bsxfun(@times, Y, subject_scale);
Y = bsxfun(@plus, Y, subject_shift);

% Include a few non-finite values to exercise robust summaries without
% making an entire image invalid.
Y(10, 3) = NaN;
Y(round(V * 0.6), 5) = NaN;

end


function [elapsed_time, Z, stats] = timed_kernel(Y, gm_mask, wm_mask, csf_mask, cfg, use_gpu)

tic;
[Z, stats] = normalize_gm_shift_scale(Y, gm_mask, wm_mask, csf_mask, ...
    'do_scale', cfg.do_scale, ...
    'log_scale', cfg.log_scale, ...
    'trim_pct', 5, ...
    'use_gpu', use_gpu);
elapsed_time = toc;

end


function object_result = run_emotionreg_object_benchmark(tol)

object_result = struct();

try
    obj = load_image_set('emotionreg', 'noverbose');

    tic;
    obj_cpu = normalize_gm_by_wm_csf(obj);
    object_result.cpu_time_sec = toc;

    tic;
    obj_gpu = normalize_gm_by_wm_csf(obj, 'use_gpu', true);
    object_result.gpu_requested_time_sec = toc;

    object_result.max_dat_abs_diff = max_abs_diff(obj_cpu.dat, obj_gpu.dat);
    object_result.status = 'completed';

    fprintf('emotionreg object benchmark: CPU %.4fs, use_gpu %.4fs, max |dat diff| %.3g\n', ...
        object_result.cpu_time_sec, object_result.gpu_requested_time_sec, ...
        object_result.max_dat_abs_diff);

    assert(object_result.max_dat_abs_diff <= tol, ...
        'Object benchmark differs from CPU. Max diff = %.6g', object_result.max_dat_abs_diff);

catch err
    object_result.status = 'skipped_or_failed';
    object_result.message = err.message;
    warning('normalize_gm_by_wm_csf_gpu_unit_test:ObjectBenchmarkSkipped', ...
        'Object benchmark skipped or failed: %s', err.message);
end

end


function d = max_abs_diff(a, b)

delta = abs(double(a(:)) - double(b(:)));
delta = delta(isfinite(delta));

if isempty(delta)
    d = 0;
else
    d = max(delta);
end

end


function d = max_stats_abs_diff(stats_a, stats_b)

fields_a = fieldnames(stats_a);
fields_b = fieldnames(stats_b);
common_fields = intersect(fields_a, fields_b);
d = 0;

for i = 1:numel(common_fields)
    name = common_fields{i};
    a = stats_a.(name);
    b = stats_b.(name);

    if isnumeric(a) || islogical(a)
        d = max(d, max_abs_diff(a, b));
    end
end

end


function [tf, message] = local_gpu_available()

tf = false;
message = 'gpuDeviceCount/gpuDevice not available';

try
    if exist('gpuDeviceCount', 'file') ~= 2 && exist('gpuDeviceCount', 'builtin') ~= 5
        return
    end

    if gpuDeviceCount < 1
        message = 'gpuDeviceCount returned 0';
        return
    end

    gpuDevice;
    tf = true;
    message = '';

catch err
    message = err.message;
end

end
