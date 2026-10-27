function [tc, meta] = hrf_extract_signature_timeseries(fmri_nii, varargin)
%HRF_EXTRACT_SIGNATURE_TIMESERIES Extract interpretable signature time-series.
% Uses apply_all_signatures on an fmri_data object.
%
% Name/value
%   'SimilarityMetric' : dotproduct (default), cosine_similarity, correlation
%   'ImageSet'         : all (default), or a named signature set accepted by apply_all_signatures
%   'SignatureName'    : optional specific signature to return; default = first available

p = inputParser;
p.addRequired('fmri_nii', @(x) ischar(x) || isstring(x) || isa(x, 'fmri_data'));
p.addParameter('SimilarityMetric', 'dotproduct', @(x) ischar(x) || isstring(x));
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

if isa(fmri_nii, 'fmri_data')
    dat = fmri_nii;
else
    dat = fmri_data(char(fmri_nii));
end
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
if istable(sig_struct)
    vn = sig_struct.Properties.VariableNames;
    tc = sig_struct.(vn{1});
elseif isstruct(sig_struct)
    f = fieldnames(sig_struct);
    tc = sig_struct.(f{1});
else
    tc = sig_struct;
end

tc = tc(:);
s = std(tc);
if s == 0 || isnan(s)
    tc = zeros(size(tc));
else
    tc = (tc - mean(tc)) ./ s;
end

meta = struct();
meta.selected_signature = selected_name;
meta.available_signatures = sig_names;
meta.similarity_metric = char(opts.SimilarityMetric);
meta.image_set = char(opts.ImageSet);
end
