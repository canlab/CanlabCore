function [tc, meta] = hrf_extract_signature_timeseries(fmri_nii, varargin)
%HRF_EXTRACT_SIGNATURE_TIMESERIES Extract interpretable signature time-series.
% Uses apply_all_signatures on an fmri_data object.
%
% Name/value
%   'SimilarityMetric' : dot_product (default), cosine_similarity, correlation
%   'ImageSet'         : all (default), or a named signature set accepted by apply_all_signatures
%   'SignatureName'    : optional specific signature to return; default = first available

p = inputParser;
p.addRequired('fmri_nii', @(x) ischar(x) || isstring(x));
p.addParameter('SimilarityMetric', 'dot_product', @(x) ischar(x) || isstring(x));
p.addParameter('ImageSet', 'all', @(x) ischar(x) || isstring(x));
p.addParameter('SignatureName', '', @(x) ischar(x) || isstring(x));
p.parse(fmri_nii, varargin{:});
opts = p.Results;

if ~exist('fmri_data', 'file')
    error('fmri_data class/function not found on path. Add CanlabCore dependencies first.');
end
if ~exist('apply_all_signatures', 'file')
    error('apply_all_signatures not found on path. Add CanlabCore dependencies first.');
end

dat = fmri_data(char(fmri_nii));
S = apply_all_signatures(dat, 'similarity_metric', char(opts.SimilarityMetric), ...
    'image_set', char(opts.ImageSet));

if ~isfield(S, 'signaturenames') || isempty(S.signaturenames)
    error('apply_all_signatures returned no signatures.');
end

sig_names = S.signaturenames;
selected_name = char(opts.SignatureName);
if isempty(selected_name)
    selected_name = sig_names{1};
elseif ~ismember(selected_name, sig_names)
    error('Requested SignatureName "%s" not in apply_all_signatures output.', selected_name);
end

sig_struct = S.(selected_name);
if isstruct(sig_struct)
    f = fieldnames(sig_struct);
    tc = sig_struct.(f{1});
else
    tc = sig_struct;
end

tc = tc(:);
tc = (tc - mean(tc)) ./ std(tc);

meta = struct();
meta.selected_signature = selected_name;
meta.available_signatures = sig_names;
meta.similarity_metric = char(opts.SimilarityMetric);
meta.image_set = char(opts.ImageSet);
end
