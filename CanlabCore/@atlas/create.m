function obj = create(obj, varargin)
% create Populate an atlas object from fieldname/value pairs.
%
% Create an object from an empty obj structure, assigning fieldname/value
% pairs as optional arguments. For known field names, the corresponding
% obj property is set to the supplied value; the 'dat' field is given
% special handling so that NaN entries are silently converted to zero.
%
% This method is used internally by the fmri_data and atlas class
% constructors.
%
% :Usage:
% ::
%
%     obj = create(obj, fieldname1, value1, fieldname2, value2, ...)
%
% :Inputs:
%
%   **obj:**
%        An atlas-class (or fmri_data) object to be populated.
%
%   **varargin:**
%        Pairs of (fieldname, value) where fieldname is a property of
%        obj. The optional keyword 'noverbose' suppresses progress
%        messages.
%
% :Outputs:
%
%   **obj:**
%        The input object with specified fields populated.
%
% :See also:
%   - atlas
%   - fmri_data

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
