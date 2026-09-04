function validate_metadata(obj)
% validate_metadata  Check internal consistency of an rsm's metadata.
%
% Throws an error if:
%   - .labels length disagrees with size(.dat, 1)
%   - .metadata_table row count disagrees with size(.dat, 1)
%   - .groupings contains indices outside 1:k
%   - .replicate_table row count disagrees with size(.dat, 3)
%   - .dat is non-square in first two dims when non-empty

if isempty(obj.dat), return; end

k = size(obj.dat, 1);

if size(obj.dat, 2) ~= k
    error('rsm:notSquare', ...
        'rsm.dat must be square in the first two dimensions; got [%d x %d].', ...
        size(obj.dat,1), size(obj.dat,2));
end

if ~isempty(obj.labels) && numel(obj.labels) ~= k
    error('rsm:labelMismatch', ...
        'numel(labels) = %d but size(dat,1) = %d.', numel(obj.labels), k);
end

if ~isempty(obj.metadata_table) && height(obj.metadata_table) ~= k
    error('rsm:metadataMismatch', ...
        'height(metadata_table) = %d but size(dat,1) = %d.', ...
        height(obj.metadata_table), k);
end

if ~isempty(obj.groupings) && isstruct(obj.groupings) && ~isempty(fieldnames(obj.groupings))
    fn = fieldnames(obj.groupings);
    for i = 1:numel(fn)
        idx = obj.groupings.(fn{i});
        if isempty(idx), continue; end
        if ~isnumeric(idx) && ~islogical(idx)
            error('rsm:groupingBadType', ...
                'groupings.%s must be numeric or logical indices.', fn{i});
        end
        if any(idx(:) < 1) || any(idx(:) > k)
            error('rsm:groupingOutOfRange', ...
                'groupings.%s contains indices outside 1:%d.', fn{i}, k);
        end
    end
end

N = size(obj.dat, 3);
if N > 1 && ~isempty(obj.replicate_table) && height(obj.replicate_table) ~= N
    error('rsm:replicateMismatch', ...
        'height(replicate_table) = %d but size(dat,3) = %d.', ...
        height(obj.replicate_table), N);
end

end
