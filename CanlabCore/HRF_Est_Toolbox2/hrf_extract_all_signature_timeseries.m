function [TC, meta] = hrf_extract_all_signature_timeseries(fmri_nii, varargin)
%HRF_EXTRACT_ALL_SIGNATURE_TIMESERIES Extract all signature-expression time series.

p = inputParser;
p.addRequired('fmri_nii', @(x) ischar(x) || isstring(x) || isa(x, 'fmri_data'));
p.addParameter('SimilarityMetric', 'dotproduct', @(x) ischar(x) || isstring(x));
p.addParameter('ImageSet', 'all', @(x) ischar(x) || isstring(x));
p.parse(fmri_nii, varargin{:});
opts = p.Results;

if ~exist('fmri_data', 'file')
    error('fmri_data class/function not found on path. Add CanlabCore dependencies first.');
end
if ~exist('apply_all_signatures', 'file')
    error('apply_all_signatures not found on path. Add CanlabCore dependencies first.');
end

if isa(fmri_nii, 'fmri_data')
    dat = fmri_nii;
else
    dat = fmri_data(char(fmri_nii));
end
S = apply_all_signatures(dat, 'similarity_metric', char(opts.SimilarityMetric), 'image_set', char(opts.ImageSet));

if ~isfield(S, 'signaturenames') || isempty(S.signaturenames)
    error('apply_all_signatures returned no signatures.');
end

names = S.signaturenames;
n_sig = numel(names);
first = local_get_signal(S.(names{1}));
n_tp = numel(first);
TC = nan(n_tp, n_sig);

for i = 1:n_sig
    v = local_get_signal(S.(names{i}));
    if numel(v) ~= n_tp
        error('Signature time series length mismatch for %s.', names{i});
    end
    TC(:, i) = local_zscore(v);
end

meta = struct('available_signatures', {names}, ...
    'similarity_metric', char(opts.SimilarityMetric), ...
    'image_set', char(opts.ImageSet));
end

function v = local_get_signal(sig_struct)
if istable(sig_struct)
    vn = sig_struct.Properties.VariableNames;
    v = sig_struct.(vn{1});
elseif isstruct(sig_struct)
    f = fieldnames(sig_struct);
    v = sig_struct.(f{1});
else
    v = sig_struct;
end
if istable(v)
    vn = v.Properties.VariableNames;
    v = v.(vn{1});
end
v = v(:);
end

function y = local_zscore(y)
y = y(:);
s = std(y);
if s == 0 || isnan(s)
    y = zeros(size(y));
else
    y = (y - mean(y)) ./ s;
end
end
