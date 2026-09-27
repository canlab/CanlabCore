function varargout = size(obj, dim)
% size  Return the size of the rsm.
%
% For an rsm array (numel(obj) > 1), returns the array shape via builtin
% (so MATLAB infrastructure that calls size() — arrayfun, disp, indexing —
% works correctly).
%
% For a scalar rsm, returns size(obj.dat) — the [k x k x N] shape of the
% similarity matrix itself.

if numel(obj) ~= 1
    % Array of rsm objects: defer to the builtin so arrayfun and the
    % MATLAB display machinery see the array shape, not the internal dat.
    if nargin < 2
        sz = builtin('size', obj);
        if nargout <= 1
            varargout{1} = sz;
        else
            for i = 1:nargout
                if i <= numel(sz), varargout{i} = sz(i);
                else, varargout{i} = 1;
                end
            end
        end
    else
        varargout{1} = builtin('size', obj, dim);
    end
    return
end

% Scalar rsm: forward to the underlying matrix
if nargin < 2
    if nargout <= 1
        varargout{1} = size(obj.dat);
    else
        sz = size(obj.dat);
        for i = 1:nargout
            if i <= numel(sz), varargout{i} = sz(i);
            else, varargout{i} = 1;
            end
        end
    end
else
    varargout{1} = size(obj.dat, dim);
end

end
