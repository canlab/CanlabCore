function obj = to_rsm(obj)
% to_rsm  Convert an RDM to an RSM (1 - dissimilarity), preserving metadata.
%
% No-op when already a similarity matrix.

if numel(obj) > 1
    obj = arrayfun(@to_rsm, obj);
    return
end

if ~obj.is_dissimilarity, return; end

obj.dat = 1 - obj.dat;
% Diagonal -> 1
k = size(obj.dat, 1);
N = size(obj.dat, 3);
for n = 1:N
    slice = obj.dat(:, :, n);
    slice(1:k+1:end) = 1;
    obj.dat(:, :, n) = slice;
end

obj.is_dissimilarity = false;
obj.history{end+1} = sprintf('%s: to_rsm (1 - dissimilarity)', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

end
