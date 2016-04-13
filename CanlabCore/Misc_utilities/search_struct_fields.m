function found_paths = search_struct_fields(search_struct, fieldname, fieldpath)
% Returns a list of all paths inside structure search_struct that match
% the fieldname (or start with it)
%
% :Usage:
% ::
%
%     found_paths = search_struct_fields(search_struct, fieldname)
%
% :Examples:
% ::
%
%    foo = [];
%    foo.foo = [];
%    foo.foo.foo = [];
%    search_struct_fields(foo, 'foo')
%    ans = 
%       'foo.foo'
%       'foo.foo.foo'
%
%    % or
%
%    search_struct_fields(SPM, 'x')
%    ans = 
%       'SPM.xX'
%       'SPM.xM'
%       'SPM.xsDes'
%       'SPM.xX.xVi'
%       'SPM.xX.xKXs'
%       'SPM.xM.xs'

    found_paths = {};
    if(~isstruct(search_struct))
        error('search_struct not a structure');
    end
    if(~exist('fieldpath', 'var') || isempty(fieldpath))
        fieldpath = inputname(1);
    end
    if(length(search_struct) > 1)
        fieldpath = [fieldpath '(:)'];
    end

    names = fieldnames(search_struct);
    matches = strmatch(fieldname, names);
    
    if(~isempty(matches))
        found_paths = cellstr([repmat([fieldpath '.'], length(matches), 1) strvcat(names{matches})]);
    end
    
    for i=1:length(names)
        if(isstruct(search_struct(1).(names{i})))
            new_fieldpath = [fieldpath '.' names{i}];
            found_paths = [found_paths; search_struct_fields(search_struct(1).(names{i}), fieldname, new_fieldpath)];
        end
    end
end
