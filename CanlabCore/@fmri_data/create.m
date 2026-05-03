function obj = create(obj, varargin)
% create Assemble an fmri_data object from fieldname/value pairs.
%
% :Usage:
% ::
%
%     obj = create(obj, 'fieldname1', value1, 'fieldname2', value2, ...)
%     obj = create(obj, ..., 'noverbose')
%
% Internal helper used by the fmri_data class constructor and other
% methods that need to populate fields of an existing object from a
% list of name/value pairs. For each name/value pair where the name
% matches a property of obj, the corresponding field is assigned.
% NaN values written to .dat are converted to 0 with a notice (unless
% 'noverbose' is supplied).
%
% :Inputs:
%
%   **obj:**
%        An existing fmri_data (or subclass) object whose fields will
%        be populated.
%
%   **varargin:**
%        Alternating fieldname / value pairs. Fieldnames not matching
%        any property of obj are silently ignored. The string
%        'noverbose' anywhere in varargin suppresses informational
%        output (e.g., the NaN-to-zero notice).
%
% :Outputs:
%
%   **obj:**
%        The input object with the requested fields populated.
%
% :See also:
%   - fmri_data (class constructor)
%
% ..
%    Programmers' notes:
%    Used in fmri_data.m class constructor.
% ..

% if 'noverbose' is entered, suppress output
verbose = isempty(strmatch('noverbose', varargin(cellfun(@ischar, varargin))));

N = fieldnames(obj);

for i = 1:length(varargin)
    if ischar(varargin{i})
        
        % Look for a field (attribute) with the input name
        wh = strmatch(varargin{i}, N, 'exact');
        
        if ~isempty(wh)
            
            obj.(varargin{i}) = varargin{i + 1};
            
            % special methods for specific fields
            switch varargin{i}
                
                case 'dat'
                    xx = isnan(obj.(varargin{i}));
                    if any(xx(:))
                        if verbose
                            fprintf('fmri_data.create: Converting %3.0f NaNs to 0s.', sum(xx(:)));
                        end
                        obj.dat(xx) = 0;
                    end
                    
            end %Switch
            
        end % ~isempty
    end
end
end % function
