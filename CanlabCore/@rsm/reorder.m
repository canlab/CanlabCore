function obj = reorder(obj, perm)
% reorder  Apply a permutation to the row/col order of an rsm.
%
% Usage
%   R2 = reorder(R, [3 1 2 ...])

if numel(obj) > 1, obj = arrayfun(@(o) reorder(o, perm), obj); return; end

k = size(obj.dat, 1);
perm = perm(:);
if numel(perm) ~= k || ~isequal(sort(perm), (1:k)')
    error('rsm:reorder:badPerm', 'perm must be a permutation of 1:%d.', k);
end

obj.dat = obj.dat(perm, perm, :);
if ~isempty(obj.labels), obj.labels = obj.labels(perm); end
if ~isempty(obj.metadata_table), obj.metadata_table = obj.metadata_table(perm, :); end

if ~isempty(obj.groupings) && ~isempty(fieldnames(obj.groupings))
    fn = fieldnames(obj.groupings);
    [~, inv_perm] = sort(perm);
    for i = 1:numel(fn)
        obj.groupings.(fn{i}) = inv_perm(obj.groupings.(fn{i}))';
    end
end

obj.history{end+1} = sprintf('%s: reorder (perm)', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

end
