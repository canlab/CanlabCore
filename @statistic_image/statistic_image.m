classdef statistic_image < image_vector
% Class definition for a statistic_image object.
% Here is a short example of some things you can do with this object:
% 
% As with any object class in this toolbox, you can create an object by
% specifying names of fields paired with values. You can also enter
% filenames when you call statistic_image and create an image by loading a
% file. Finally, some methods performed on other objects (e.g., predict,
% regress, and ttest for fmri_data objects) will return statistic_image
% objects.
% For example, the line of code below creates a map of random t-values:
% t = statistic_image('dat', rand(1000, 1), 'type', 't', 'dfe', 29);
%
% If you want to use the full functionality of the object type without
% errors, however, a .volInfo field will need to be attached with
% information about the image space.  The easiest way to create this is to
% load an existing image, which will automatically read in the space along 
% with other info:
%
% Start with the name of a statistic image we're interested in, with arbitrary values
% name = 'salientmap.nii';   
% 
% % Load it into a statistic_image object
% img = statistic_image('image_names', name);
% 
% % Plot a histogram of its values
% figure; histogram(img)
% 
% % Threshold the image values; save values < 0 or > 5
% img = threshold(img, [0 5], 'raw-outside');
% 
% % Plot a montage and spm-orthviews of the thresholded image
% montage(img);
% orthviews(img);
% 
% % Remove the threshold; include all values between -Inf and Inf
% img = threshold(img, [-Inf Inf], 'raw-between');
% 
% % Show the orthviews plot of the unthresholded image
% orthviews(img);
%
% You can force an image to be of a particular type, e.g., 'T', in which
% case p-values will automatically be added. This is useful for
% thresholding. 
% e.g., for an SPM T-map: 
% wh = 1; 
% load SPM
% img = sprintf('spmT_00%02d.img', wh(1))
% t = statistic_image('image_names', img, 'type', 'T', 'dfe', SPM.xX.erdf);
%
% If 'type' is 'robreg', then assumes user is in a robust regression
% directory and will load in the robust_beta_000X.img and load in
% robust_p_000X.img as the p values.  X is assumed = 1, but can be
% overloaded by passing in a dat_descrip field. i.e.
% deltadon=statistic_image('dat_descrip', 4, 'type', 'robreg')


    properties
        % properties of parent class are also inherited
        
        type
        
        p
        
        p_type
        
        ste
        
        threshold
        
        thr_type
        
        sig
        
        N
        
        dfe
        
    end
    
    methods
        
        function obj = statistic_image(varargin)
            
            % Enter fieldname', value pairs in any order to create class
            % instance
            
            % ---------------------------------
            % Create empty image_vector object, and return if no additional
            % arguments
            % ---------------------------------
            obj.type = 'generic';
            %obj.X = []; alread defined
            
            obj.p = [];
            obj.p_type = [];
            obj.ste = [];
            obj.threshold = [];
            obj.thr_type = [];
            obj.sig = [];
            obj.dfe = [];
            obj.N = [];
            
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
            
            % if loading in a robust regression image, overload image_names
            obj.image_names = sprintf('rob_beta_%04d.img',obj.dat_descrip);  
            
            % load data, if we don't have data yet
            % But we do have valid image names.
            % ----------------------------------------
            
            obj = check_image_filenames(obj);
            
            if isempty(obj.dat) && any(obj.files_exist)
                obj = read_from_file(obj);

            elseif isempty(obj.dat)
                disp('Warning: .dat is empty and files cannot be found.  No image data in object.');
            end
            
            
            % 6/22/13 Tor Added to enforce consistency in objects across usage cases
            if isempty(obj.sig) || (numel(obj.sig) == 1 && ~obj.sig)
                obj.sig = true(size(obj.dat));
            end

            % type - handle specific forced types
            
            switch obj.type    
                case {'T', 't'}
                    obj.type = 'T';
                    
                    if isempty(obj.dfe)
                        error('If forcing object type = ''T'', enter dfe in call to object constructor.')
                    end
                    
                    obj.p = 2 * (1 - tcdf(abs(obj.dat), obj.dfe));
                    obj.p_type = '2-tailed P-value from input dfe';
                    
                case 'p'
                    obj.type = 'p';
                    obj.p = obj.dat;
                
                case 'robreg'
                    obj.type = 'robust regression:  .dat is beta values, .p is robust p values';
                    
                    pimg = image_vector('image_names', sprintf('rob_p_%04d.img',obj.dat_descrip));
                    obj.p = pimg.dat;
                    obj.p_type = 'robust';
                    obj.dat_descrip = ['robust regression for covariate number ' num2str(obj.dat_descrip)];
                    
                case 'generic'
                    % do nothing
                    
                otherwise
                    error('User forced unknown statistic image type.')
            end
            
            
        end % class constructor
        
    end  % methods
    
end  % classdef

