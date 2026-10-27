function obj = fisher_z(obj)
% fisher_z  Apply atanh to .dat, in-place value transform.
%
% Only sensible for correlation-style RSMs (values in [-1, 1]). For metrics
% outside that range, returns unchanged with a warning.

if numel(obj) > 1, obj = arrayfun(@fisher_z, obj); return; end

if obj.is_dissimilarity || ~ismember(lower(obj.metric), {'correlation','spearman','cosine'})
    warning('rsm:fisher_z:unusualMetric', ...
        'fisher_z applied to a non-correlation rsm (metric=%s); proceeding but check values.', obj.metric);
end

% Clip to avoid Inf at exactly ±1
d = obj.dat;
d(d >  0.9999999) =  0.9999999;
d(d < -0.9999999) = -0.9999999;
obj.dat = atanh(d);

obj.history{end+1} = sprintf('%s: fisher_z (atanh)', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

end
