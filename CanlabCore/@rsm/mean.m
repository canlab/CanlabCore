function obj = mean(obj, varargin)
% mean  Average an rsm across the replicate axis (3rd dim).
%
% Usage
%   R_avg = mean(R)                    % nanmean across dim 3
%   R_avg = mean(R, 'fisherz', true)   % Fisher-z transform before averaging, untransform after
%   R_avg = mean(R, 'omitnan', false)  % use mean() instead of nanmean()
%
% Returns a new rsm with size [k x k x 1] and a single-row replicate_table
% (or empty) indicating the aggregation.

if numel(obj) > 1, error('rsm:mean:nonScalar', 'mean() expects a scalar rsm; use arrayfun for arrays.'); end

p = inputParser;
p.addParameter('fisherz', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('omitnan', true,  @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

d = obj.dat;
if size(d, 3) <= 1, return; end   % nothing to average

use_fz = logical(p.Results.fisherz) && ~obj.is_dissimilarity && ...
         ismember(lower(obj.metric), {'correlation','spearman','cosine'});

if use_fz
    d_clip = d;
    d_clip(d_clip >  0.9999999) =  0.9999999;
    d_clip(d_clip < -0.9999999) = -0.9999999;
    d = atanh(d_clip);
end

if p.Results.omitnan
    m = mean(d, 3, 'omitnan');
else
    m = mean(d, 3);
end

if use_fz
    m = tanh(m);
end

obj.dat = m;

if ~isempty(obj.replicate_table)
    % Collapse the replicate table to a single row describing the aggregation
    obj.replicate_table = table({'mean'}, 'VariableNames', {'aggregation'});
end

obj.level = 'collapsed';
obj.history{end+1} = sprintf('%s: mean across replicate axis (fisherz=%d)', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), use_fz);

end
