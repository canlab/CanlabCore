function tf = isempty(obj)
% isempty  True if obj has no elements OR (scalar case) obj.dat is empty.

n = builtin('numel', obj);

if n == 0
    tf = true;
    return
end

if n > 1
    tf = false(builtin('size', obj));
    for i = 1:n
        tf(i) = isempty(obj(i).dat);
    end
    return
end

tf = isempty(obj.dat);

end
