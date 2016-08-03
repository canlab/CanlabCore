function obj = surface(obj, varargin)
% Adds surfaces of brain to figure
%
% :Usage:
% ::
%
%     obj = surface(obj, varargin)
%
% :Inputs:
%
%   **obj:**
%        fmridisplay
%
% :Outputs:
%
%   **obj:**
%        an fmridisplay object
%
% :Properties:
%
%  - overlay: ''
%  - SPACE: ''
%  - activation_maps: {}
%  - montage: {}
%  - surface: {[1x1 struct]}
%  - orthviews: {}
%  - history: {}
%  - history_descrip: []
%  - additional_info: ''
%
% :Examples:
%
%     o2 = surface(o2, axes, [0.15 0.28 .15 1], 'direction', 'hires right', 'orientation', 'lateral');
%
% See help fmridisplay

if nargin == 0 || isempty(obj) || ~isa(obj, 'fmridisplay')
    obj = fmridisplay;
end
% initialize, if nothing passed in; but you would have to call overloaded
% method, fmridisplay.montage, to invoke this.


if isempty(obj.SPACE) || ~isstruct(obj.SPACE.V)
    error('fmridisplay is not initialized correctly. run obj = fmridisplay; first and then pass in to this method.')
    
end

ax = gca;
dir = 'hires right';
orn = 'lateral';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            
            case {'direction'}, dir = varargin{i + 1};
            
            case {'orientation'}, orn = varargin{i + 1};

            case {'axes'}, ax = varargin{i + 1};

                %otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

axh = axes('Position', ax);
% axh = axes(ax);
if strcmp(dir, 'hires left')
    h = addbrain('hires left');
    if strcmp(orn, 'medial')
        view(270, 0);
    else
        view(90, 0);
    end
elseif strcmp(dir, 'left')
    h = addbrain('left');
    if strcmp(orn, 'medial')
        view(270, 0);
    else
        view(90, 0);
    end
elseif strcmp(dir, 'hires right')
    h = addbrain('hires right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
elseif strcmp(dir, 'right')
    h = addbrain('right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
elseif strcmp(dir, 'flat right')
    h = addbrain('flat right');
    
elseif strcmp(dir, 'flat right')
    h = addbrain('right');
    
elseif strcmp(dir, 'surface left')
    h = addbrain('surface left');
        
    if strcmp(orn, 'medial')
        view(270, 0);
    else
        view(90, 0);
    end
    
elseif strcmp(dir, 'surface right')
    h = addbrain('surface right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
   
 
end

if strcmp(h(1).FaceColor,'interp')
    set(h,'FaceAlpha', 1);
else
    set(h, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
end

lightRestoreSingle; 
axis image;
axis off; 
lighting gouraud;
material dull;


% register info in fmridisplay object
% 'axis_handles', axh, 
obj.surface{end + 1} = struct('axis_handles', axh, 'direction', dir, 'orientation', orn, 'object_handle', h);

end % main function

