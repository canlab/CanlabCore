function obj = subset(obj, idx)
% subset  Restrict an rsm to a subset of conditions.
%
% Usage
%   R2 = subset(R, 1:8)
%   R2 = subset(R, [true false true ...])
%   R2 = subset(R, 'hot')               % grouping name (uses .groupings.hot)
%
% Preserves replicate axis; subsets labels, metadata_table, and groupings
% accordingly.

if numel(obj) > 1, obj = arrayfun(@(o) subset(o, idx), obj); return; end

% Resolve grouping name to indices
if ischar(idx) || isstring(idx)
    name = char(idx);
    if ~isfield(obj.groupings, name)
        error('rsm:subset:unknownGrouping', 'No grouping named "%s".', name);
    end
    idx = obj.groupings.(name);
end

if islogical(idx), idx = find(idx); end
idx = idx(:);
k = size(obj.dat, 1);
if any(idx < 1) || any(idx > k)
    error('rsm:subset:outOfRange', 'subset indices must be in 1:%d.', k);
end

obj.dat = obj.dat(idx, idx, :);

if ~isempty(obj.labels), obj.labels = obj.labels(idx); end
if ~isempty(obj.metadata_table), obj.metadata_table = obj.metadata_table(idx, :); end

% Reindex groupings (drop any whose indices are not all in the subset)
if ~isempty(obj.groupings) && ~isempty(fieldnames(obj.groupings))
    fn = fieldnames(obj.groupings);
    new_groupings = struct();
    for i = 1:numel(fn)
        old_idx = obj.groupings.(fn{i});
        [tf, loc] = ismember(old_idx, idx);
        if all(tf)
            new_groupings.(fn{i}) = loc(:)';
        end  % otherwise drop this grouping
    end
    obj.groupings = new_groupings;
end

obj.history{end+1} = sprintf('%s: subset (k=%d -> %d)', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), k, numel(idx));

end
