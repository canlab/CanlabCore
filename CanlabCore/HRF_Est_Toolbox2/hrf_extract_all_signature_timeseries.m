function [TC, meta] = hrf_extract_all_signature_timeseries(fmri_nii, varargin)
%HRF_EXTRACT_ALL_SIGNATURE_TIMESERIES Extract all signature-expression time series.

p = inputParser;
p.addRequired('fmri_nii', @(x) ischar(x) || isstring(x));
p.addParameter('SimilarityMetric', 'dot_product', @(x) ischar(x) || isstring(x));
p.addParameter('ImageSet', 'all', @(x) ischar(x) || isstring(x));
p.parse(fmri_nii, varargin{:});
opts = p.Results;

if ~exist('fmri_data', 'file')
    error('fmri_data class/function not found on path. Add CanlabCore dependencies first.');
end
if ~exist('apply_all_signatures', 'file')
    error('apply_all_signatures not found on path. Add CanlabCore dependencies first.');
end

dat = fmri_data(char(fmri_nii));
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
    v = v(:);
    TC(:, i) = (v - mean(v)) ./ std(v);
end

meta = struct('available_signatures', {names}, ...
    'similarity_metric', char(opts.SimilarityMetric), ...
    'image_set', char(opts.ImageSet));
end

function v = local_get_signal(sig_struct)
if isstruct(sig_struct)
    f = fieldnames(sig_struct);
    v = sig_struct.(f{1});
else
    v = sig_struct;
end
v = v(:);
end
