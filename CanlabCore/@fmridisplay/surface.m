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
    
elseif strcmp(dir, 'flat left')
    h = addbrain('flat left');
    
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

    
elseif strcmp(dir, 'inflated right')
    h = addbrain('inflated right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'inflated left')
    h = addbrain('inflated left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
   
elseif strcmp(dir, 'hcp inflated right')
    h = addbrain('hcp inflated right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'hcp inflated left')
    h = addbrain('hcp inflated left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
   
elseif strcmp(dir, 'hcp sphere right')
    h = addbrain('hcp sphere right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'hcp sphere left')
    h = addbrain('hcp sphere left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end

elseif strcmp(dir, 'freesurfer sphere right')
    h = addbrain('freesurfer sphere right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'freesurfer sphere left')
    h = addbrain('freesurfer sphere left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end

elseif strcmp(dir, 'freesurfer white right')
    h = addbrain('freesurfer white right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'freesurfer white left')
    h = addbrain('freesurfer white left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end

elseif strcmp(dir, 'freesurfer inflated right')
    h = addbrain('freesurfer inflated right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'freesurfer inflated left')
    h = addbrain('freesurfer inflated left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end

elseif strcmp(dir, 'MNI152NLin2009cAsym white right')
    % for internal development use
    h = addbrain('MNI152NLin2009cAsym white right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'MNI152NLin2009cAsym white left')
    % for internal development use
    h = addbrain('MNI152NLin2009cAsym white left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end

elseif strcmp(dir, 'MNI152NLin2009cAsym midthickness right')
    % for internal development use
    h = addbrain('MNI152NLin2009cAsym midthickness right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'MNI152NLin2009cAsym midthickness left')
    % for internal development use
    h = addbrain('MNI152NLin2009cAsym midthickness left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end

elseif strcmp(dir, 'MNI152NLin2009cAsym pial right')
    % for internal development use
    h = addbrain('MNI152NLin2009cAsym pial right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'MNI152NLin2009cAsym pial left')
    % for internal development use
    h = addbrain('MNI152NLin2009cAsym pial left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end

elseif strcmp(dir, 'MNI152NLin6Asym white right')
    % for internal development use
    h = addbrain('MNI152NLin6Asym white right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'MNI152NLin6Asym white left')
    % for internal development use
    h = addbrain('MNI152NLin6Asym white left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end

elseif strcmp(dir, 'MNI152NLin6Asym midthickness right')
    % for internal development use
    h = addbrain('MNI152NLin6Asym midthickness right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'MNI152NLin6Asym midthickness left')
    % for internal development use
    h = addbrain('MNI152NLin6Asym midthickness left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end

elseif strcmp(dir, 'MNI152NLin6Asym pial right')
    % for internal development use
    h = addbrain('MNI152NLin6Asym pial right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'MNI152NLin6Asym pial left')
    % for internal development use
    h = addbrain('MNI152NLin6Asym pial left');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end

elseif strcmp(dir, 'MNI152NLin6Asym sphere right')
    % for internal development use
    h = addbrain('MNI152NLin6Asym sphere right');
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'MNI152NLin6Asym sphere left')
    % for internal development use
    h = addbrain('MNI152NLin6Asym sphere left');
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

