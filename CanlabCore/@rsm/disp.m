function disp(obj)
% disp  Custom textual display of an rsm in the command window.

n = builtin('numel', obj);

if n == 0
    fprintf('  empty rsm array (0 elements)\n\n');
    return
end

if n > 1
    sz = builtin('size', obj);
    fprintf('  [%s] array of rsm objects\n\n', ...
        strjoin(arrayfun(@num2str, sz, 'UniformOutput', false), 'x'));
    max_show = min(n, 10);
    for i = 1:max_show
        fprintf('  (%d) ', i); display_one(obj(i));
    end
    if n > max_show
        fprintf('  ... and %d more (use obj(i) to inspect, or get_by_label / select)\n\n', n - max_show);
    end
    return
end

display_one(obj);

end


function display_one(obj)

if isempty(obj.dat)
    fprintf('  empty rsm  (metric=%s)\n', obj.metric);
    return
end

k = size(obj.dat, 1);
N = size(obj.dat, 3);

kind = 'RSM'; if obj.is_dissimilarity, kind = 'RDM'; end

fprintf('  %s [%d x %d]', kind, k, k);
if N > 1
    fprintf(' x %d (%s)', N, obj.level);
end
fprintf('  metric=%s', obj.metric);
if ~strcmp(obj.whitened.level, 'none')
    fprintf('  whitened=%s', obj.whitened.level);
end
fprintf('\n');

if ~isempty(obj.labels)
    if numel(obj.labels) > 6
        show = strjoin([obj.labels(1:3)' {'...'} obj.labels(end-1:end)'], ', ');
    else
        show = strjoin(obj.labels, ', ');
    end
    fprintf('    labels:    %s\n', show);
end

if ~isempty(obj.groupings) && ~isempty(fieldnames(obj.groupings))
    fprintf('    groupings: %s\n', strjoin(fieldnames(obj.groupings), ', '));
end

fprintf('    source:    %s\n', obj.source);
fprintf('\n');

end
