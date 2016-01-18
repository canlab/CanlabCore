function a = append(a, b, field)
% Appends second parameter to the end of the first. If field is set,
% writes into the field instead of directly into the variable.
%
% :Usage:
% ::
%
%     a = append(a, data, [field])
%
% This function solely exists because Matlab has no ability to create
% empty objects (e.g. with one dimension being 0) and thus, trying to append
% to an empty array results in type cast errors, since empty matrices are
% double by default.

    if(exist('field', 'var'))
        if(isempty(a))
            a(1).(field) = b;
        else
            a(end+1).(field) = b;
        end
    else
        if(isempty(a))
            a = b;
        else
            a(end+1) = b;
        end
    end
end
