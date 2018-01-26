function obj = create(obj, varargin)
% Create an object from an empty obj structure, assigning fieldname/value
% pairs as optional arguments.
%
% :Usage:
% ::
%
%     [obj = create(obj, varargin)
%
% Used in fmri_data.m class constructor.

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
