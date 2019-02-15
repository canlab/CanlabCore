% statistic_image: An object that allows for storage and manipulation of
% images containing t-values, p-values, and thresholded results
% - multiple images can be stored in a single object
% - this is a subclass of image_vector and inherits all its methods
% - easy thresholding without changing underlying image (threshold method)
% - easy visualization (orthviews, surface, and montage methods)
%
% Usage:
% -------------------------------------------------------------------------
%  object constructor: this function creates an object and populates it
%  with data.  Entering image name(s) loads in the data.
% 
%  obj = statistic_image(varargin)
%  obj = statistic_image(nifti/img filename, or allowed fields followed by values)
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% Author and copyright information:
% -------------------------------------------------------------------------
%     Copyright (C) 2010 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Inputs:
% -------------------------------------------------------------------------
% :properties: any valid attribute/property of statistic_image, followed by
% data values to store
%   - see properties(statistic_image)
%
% :image names: a character array with one or more filenames for
% Analyze/NIFTI images (.img/nii).  
%
%
% Outputs:
% -------------------------------------------------------------------------
% obj           A statistic_image object
%
% Examples:
% -------------------------------------------------------------------------
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
% Read in a p-value image as a p-value object type:
% p = which('nonnoc_v11_4_137subjmap_weighted_pvalue.nii');
% pstat = statistic_image(p, 'type', 'p');
% pstat = threshold(pstat, .05, 'fdr');
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
%
%
% See also:
% image_vector.m
% fmri_data.m
% region.m
% fmridisplay.m
%
classdef statistic_image < image_vector

    properties

        type        % String with image type: 'generic', 't', 'p', 'robreg'
        p           % Matrix of p-values for images, [voxels x images]
        p_type      % String with source info for p-values
        ste         % Matrix of standard error values for images, [voxels x images]
        threshold   % Latest statistical threshold applied
        thr_type    % Information about statistical threshold applied
        sig         % Logical matrix of which voxels are significant, [voxels x images]
        N           % Sample size
        dfe         % Error degrees of freedom for test
        
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
            
            doverbose = false;
            
            wh = strcmp(varargin, 'verbose');
            if any(wh), doverbose = true; varargin(wh) = []; end
            
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
                        % Assume is char array of image to load
                        % could also use obj = read_from_file(obj);
                        
                        if doverbose
                            fprintf('Reading from file: %s\n', varargin{i});
                            tmp = fmri_data(varargin{i});
                        else
                            tmp = fmri_data(varargin{i}, 'noverbose');
                        end
                        
                        obj.dat = tmp.dat;
                        obj.volInfo = tmp.volInfo;
                        obj.image_names = tmp.image_names;
                        obj.fullpath = tmp.fullpath;
                        obj.files_exist = tmp.files_exist;
                        obj.history = tmp.history;
                        %warning('inputargs:BadInput', sprintf('Unknown field: %s', varargin{i}));
                        
                    end
                end % string input
            end % process inputs
            
            % load data, if we don't have data yet
            % But we do have valid image names.
            % ----------------------------------------
            
            % if loading in a robust regression image, overload image_names
            %    Yoni originally added this part (2014-07-24), and Wani
            %    added if-statement to fix a bug related to this part (2014-08-26). 
            if strcmp(obj.type, 'robreg')
                if isempty(obj.dat_descrip), obj.dat_descrip=1; end %by default, assume user wants rob img 1 if none entered
                obj.image_names = sprintf('rob_tmap_%04d.img',obj.dat_descrip); 
            end
            
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
                
                case 'r'
                    obj.type = 'r';
                    
                case 'robreg'
                    obj.type = 'robust regression:  .dat is beta values, .p is robust p values';
                    
                    pimg = image_vector('image_names', sprintf('rob_p_%04d.img',obj.dat_descrip));
                    bimg = image_vector('image_names', sprintf('rob_beta_%04d.img',obj.dat_descrip));
                    timg = image_vector('image_names', sprintf('rob_tmap_%04d.img',obj.dat_descrip));
                    
                    
                    obj.p = pimg.dat;
                    obj.p_type = 'robust';
                    obj.dat = timg.dat;%bimg.dat;
                    obj.dat_descrip = ['robust regression for covariate number ' num2str(obj.dat_descrip)];            
                    
                case 'generic'
                    % do nothing
                    
                otherwise
                    error('User forced unknown statistic image type.')
            end
            
            
        end % class constructor
        
    end  % methods
    
end  % classdef

