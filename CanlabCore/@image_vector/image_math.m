function obj_out = image_math(obj1, varargin)
% Perform simple mathematical and boolean operations on image objects
%
% :Usage:
% ::
%
%    obj_out = image_math(obj1, [optional inputs, e.g., a 2nd object, keywords])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% :Inputs:
%
%   **obj1:**
%        An image_vector object
%
% :Optional Inputs:
%
%   **obj2:**
%        An additional image_vector object
%   **{'add', 'plus'}:**
%        Keyword to perform image-wise addition of images in obj1
%                              and obj2.  Assumes these are paired/matched objects.
%   **{'subtract', 'minus'}:**
%        Keyword to perform image-wise subtraction of images
%                              in obj1 and obj2
%   **{'cat', 'concatenate'}:**
%        Concatenate obj1 and obj2 image-wise.  Requires same
%                              number of voxels in both image sets.  Returns effects
%                              codes of 1, -1 in obj_out.Y.
%   **{'power'}:**
%        Keyword to raise data to power element-wise; obj.dat = obj.dat.^b;
%                              Followed by exponent to apply (b)
%
% :Outputs:
%
%   **obj_out:**
%        The result - an image_vector object
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2015  Tor Wager
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
% ..

% ..
% DEFAULTS AND INPUTS
% ..

keyword = '';  % initalize optional variables to default values here.
% keyword. should enter one...
obj2 = [];
my_exponent = [];

% optional inputs with default values
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case {'add', 'plus'}, keyword = 'plus'; varargin{i} = [];
            case {'subtract', 'minus'}, keyword = 'minus'; varargin{i} = [];
            case {'cat', 'concatenate'}, keyword = 'cat'; varargin{i} = [];
                
            case 'power', keyword = 'power'; varargin{i+1} = my_exponent;
                %case 'basistype', basistype = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
        
    elseif isa(varargin{i}, 'image_vector')
        
        obj2 = varargin{i}; varargin{i} = [];
        
    end
end

% -------------------------------------------------------------------------
% DATA CHECKS
% -------------------------------------------------------------------------

switch keyword
    
    case ''
        disp('Nothing to do.')
        return
        
    case 'power'
        %Check if image_vector object
        if ~isa(obj1,'image_vector') 
            error('Input Data is not an image_vector object')
        end 
        
    case {'plus', 'minus', 'cat'}
        n1 = size(obj1.dat, 2);
        n2 = size(obj2.dat, 2);
        
        y1 = ones(n1, 1);
        y2 = ones(n2, 1);
        
        %Check if image_vector object
        if ~isa(obj1,'image_vector') || ~isa(obj2,'image_vector')
            error('Input Data is not an image_vector object')
        end
        
        %Check number of rows
        if size(obj1.dat,1)~=size(obj2.dat,1)
            error('number of voxels is different between objects.')
        end
        
        if strcmp(keyword, 'cat')
            % we are done
    
        else
            
            % The rest is for plus/minus, which assume paired images
            
            %Check number of columns
            if n1 ~= n2
                error('Sizes of objects do not match.');
            end
            
        end
        
    otherwise
        warning(['Unknown keyword:' keyword]);
        return
        
end % keyword ; data checks


% -------------------------------------------------------------------------
% RUN OPERATION
% -------------------------------------------------------------------------

switch keyword
    
    case ''
        disp('Nothing to do.')
        return
        
    case 'plus'
        
        obj_out = obj1;
        obj_out.dat = obj_out.dat + obj2.dat;
        obj_out.history{end+1} = 'Image-wise addition operation by image_math';
        obj_out.dat_descrip = cell(1, 3);
        obj_out.dat_descrip{1} = 'Names of images added in next cells, 1st set plus 2nd';
        obj_out.dat_descrip{2} = obj1.fullpath;
        obj_out.dat_descrip{3} = obj2.fullpath;
        obj_out.image_names = [];
        obj_out.fullpath = [];
        obj_out.files_exist = false;
        
    case 'minus'
        
        obj_out = obj1;
        obj_out.dat = obj_out.dat - obj2.dat;
        obj_out.history{end+1} = 'Image-wise subtraction operation by image_math';
        obj_out.dat_descrip = cell(1, 3);
        obj_out.dat_descrip{1} = 'Names of images subtracted in next cells, 1st set minus 2nd';
        obj_out.dat_descrip{2} = obj1.fullpath;
        obj_out.dat_descrip{3} = obj2.fullpath;
        obj_out.image_names = [];
        obj_out.fullpath = [];
        obj_out.files_exist = false;
        
    case 'cat'
        
        obj_out = obj1;
        obj_out.dat = [obj_out.dat obj2.dat];
        obj_out.history{end+1} = 'Image-wise concatenation operation by image_math';
        obj_out.image_names = char(obj1.image_names, obj2.image_names);
        obj_out.fullpath = char(obj1.fullpath, obj2.fullpath);
        if isfield(obj_out, 'Y')
            obj_out.Y = [y1; -y2];
            obj_out.Y_descrip = 'Effects codes for image set A (1) and B (-1) added by image_math';
        end
    case 'power'
        obj_out = obj1;
        obj_out.dat = obj_out.dat .^ my_exponent;
        obj_out.history{end+1} = sprintf('Raised to %3.0f element-wise by image_math', my_exponent);

    otherwise
        warning(['Unknown keyword:' keyword]);
        return
        
end % switch keyword; operation


end % main function


