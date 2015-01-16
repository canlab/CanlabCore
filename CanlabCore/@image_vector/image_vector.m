classdef image_vector
    
    properties
        
        dat % generic name for data
        dat_descrip 
        volInfo  
           
        removed_voxels = logical(0);
        removed_images = logical(0);
                
        image_names    
        fullpath  
        files_exist
        history
         
    end
    
    methods
        
        % Class constructor
        function obj = image_vector(varargin)
            % Enter fieldname', value pairs in any order to create class
            % instance
            
            % ---------------------------------
            % Create empty image_vector object, and return if no additional
            % arguments
            % ---------------------------------
            
            obj.dat = [];
            obj.volInfo = [];
            obj.image_names = [];
            obj.fullpath = char([]);
            obj.files_exist = false;
            obj.history = {};
            
            % The code below can be generic to any class definition
            % It parses 'fieldname', value pairs of inputs
            % and returns a warning if unexpected strings are found.
            
            if nargin == 0
                return
            end
            
            % all valid fieldnames
            valid_names = fieldnames(obj);
            
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    
                    % Look for a field (attribute) with the input name
                    wh = strmatch(varargin{i}, valid_names, 'exact');
                    
                    % behaviors for valid fields
                    if ~isempty(wh)
                        
                        obj.(varargin{i}) = varargin{i + 1};
                        
                        % eliminate strings to prevent warnings on char
                        % inputs
                        if ischar(varargin{i + 1})
                            varargin{i + 1} = [];
                        end
                        
                        % special methods for specific fields
                        switch varargin{i}
                            
                        end
                        
                    else
                        warning('inputargs:BadInput', sprintf('Unknown field: %s', varargin{i}));
                        
                    end
                end % string input
            end % process inputs
            
            % load data, if we don't have data yet
            % But we do have valid image names.
            % ----------------------------------------
            
            obj = check_image_filenames(obj);
            
            if isempty(obj.dat) && any(obj.files_exist)
                obj = read_from_file(obj);
                
            elseif isempty(obj.dat)
                disp('Warning: .dat is empty and files cannot be found.  No image data in object.');
            end
            
                
        end % class constructor
        
    end  % methods
    
end  % classdef

