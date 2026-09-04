function obj = to_rdm(obj)
% to_rdm  Convert an RSM to an RDM (1 - similarity), preserving metadata.
%
% No-op when already a dissimilarity matrix.

if numel(obj) > 1
    obj = arrayfun(@to_rdm, obj);
    return
end

if obj.is_dissimilarity, return; end

obj.dat = 1 - obj.dat;
% Diagonal -> 0
k = size(obj.dat, 1);
N = size(obj.dat, 3);
for n = 1:N
    slice = obj.dat(:, :, n);
    slice(1:k+1:end) = 0;
    obj.dat(:, :, n) = slice;
end

obj.is_dissimilarity = true;
obj.history{end+1} = sprintf('%s: to_rdm (1 - similarity)', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

end
