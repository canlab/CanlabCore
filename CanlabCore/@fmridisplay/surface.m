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

% Multi-surface keywords (e.g. 'foursurfaces', 'foursurfaces_hcp') add a SET
% of surface views to THIS SAME object, laid out in the current figure. Each
% becomes its own registered view, so blobs/refresh act on all of them
% together. (fmridisplay is a handle class, so han is not replaced.)
multi = '';
for kw = {'foursurfaces_hcp', 'foursurfaces_freesurfer', 'foursurfaces'}
    if any(strcmp(varargin, kw{1})), multi = kw{1}; break; end
end
if ~isempty(multi)
    obj = add_multiple_surfaces(obj, multi, varargin{:});
    return
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

if isa(ax,'matlab.graphics.axis.Axes')
    axh = ax;
    axes(axh);
else
    axh = axes('Position', ax);
end
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
    h = addbrain('hcp inflated right',varargin{:});
    if strcmp(orn, 'medial')
        view(90, 0);
    else
        view(270, 0);
    end
    
elseif strcmp(dir, 'hcp inflated left')
    h = addbrain('hcp inflated left',varargin{:});
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

elseif strcmp(dir, 'bigbrain left')
    h = addbrain('bigbrain');
    view(-137, 18); lightRestoreSingle;

elseif strcmp(dir, 'bigbrain right')
    h = addbrain('bigbrain');
    view(137, 18); lightRestoreSingle;

elseif strcmp(dir, 'left_cutaway')
    h = addbrain('left_cutaway');

elseif strcmp(dir, 'right_cutaway')
    h = addbrain('right_cutaway');

elseif strcmp(dir, 'brainstem left')
    h = addbrain('midbrain_group');
    h = [h addbrain('rvm')];
    h = [h addbrain('lc')];
    h = [h addbrain('brainstem')];
    h = [h addbrain('thalamus_group')];
    h = [h addbrain('pbn')];
    h = [h addbrain('rn')];
    h = [h addbrain('pag')];
    % h = [surface_handles addbrain('caudate')];
    view(-137, 18); lightRestoreSingle;

elseif strcmp(dir, 'brainstem right')
    h = addbrain('midbrain_group');
    h = [h addbrain('rvm')];
    h = [h addbrain('lc')];
    h = [h addbrain('brainstem')];
    h = [h addbrain('thalamus_group')];
    h = [h addbrain('pbn')];
    h = [h addbrain('rn')];
    h = [h addbrain('pag')];
    % h = [surface_handles addbrain('caudate')];
    view(137, 18); lightRestoreSingle;

elseif strcmp(dir, 'caudate left')
    h = addbrain('caudate');
    h = [h addbrain('put')];
    h = [h addbrain('gp')];
    h = [h addbrain('nacc')];
    h = [h addbrain('sn')];
    view(-137, 18); lightRestoreSingle;

elseif strcmp(dir, 'caudate right')
    h = addbrain('caudate');
    h = [h addbrain('put')];
    h = [h addbrain('gp')];
    h = [h addbrain('nacc')];
    h = [h addbrain('sn')];
    view(137, 18); lightRestoreSingle;

else
    h = addbrain(dir);

end





if strcmp(h(1).FaceColor,'interp')
    set(h,'FaceAlpha', 1);
else
    set(h, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
end
%[
lightRestoreSingle; 
axis image;
axis off; 
lighting gouraud;
material dull;


% register info in fmridisplay object
% 'axis_handles', axh,
obj.surface{end + 1} = struct('axis_handles', axh, 'direction', dir, 'orientation', orn, 'object_handle', h);

% Pull-in: if blob layers already exist on this object, render them onto the
% surface we just added, so a surface added AFTER addblobs still shows the
% blobs and add/remove/refresh act on all views together (design notes 4.2).
new_surf_idx = numel(obj.surface);
for k = 1:numel(obj.activation_maps)
    obj = render_layer_surfaces(obj, k, new_surf_idx);
end

end % main function


function obj = add_multiple_surfaces(obj, multi, varargin)
% Add a canonical set of four surface views (L/R lateral + L/R medial) to the
% same fmridisplay object, arranged 2x2 in a dedicated figure. Each is added
% via a normal single-surface surface() call (handle semantics keep obj the
% same), which also triggers blob pull-in per surface.

% Ensure a dedicated figure for the 2x2 layout. If the current figure already
% has content (e.g. a montage), the four surfaces would otherwise be drawn over
% it (and appear to render nothing). Reuse an empty current figure (e.g. the
% user typed `figure;` first); otherwise open a new one.
curfig = get(groot, 'CurrentFigure');
if isempty(curfig) || ~isempty(findobj(curfig, 'Type', 'axes'))
    figure;
end

rest = varargin;
rest(strcmp(rest, multi)) = [];                       % drop the multi keyword
for kw = {'direction', 'orientation', 'axes'}         % drop any single-view controls
    wh = find(strcmp(rest, kw{1}));
    for j = sort(wh, 'descend'), rest(j:j+1) = []; end
end

switch multi
    case 'foursurfaces_hcp'
        Lsurf = 'hcp inflated left';        Rsurf = 'hcp inflated right';
    case 'foursurfaces_freesurfer'
        Lsurf = 'freesurfer inflated left'; Rsurf = 'freesurfer inflated right';
    otherwise % 'foursurfaces'
        Lsurf = 'hires left';               Rsurf = 'hires right';
end

% {direction, orientation, axes-position} for a 2x2 layout:
% top row = lateral views, bottom row = medial views; left col = left hemi.
configs = { ...
    {Lsurf, 'lateral', [0.03 0.50 0.45 0.45]}, ...
    {Rsurf, 'lateral', [0.52 0.50 0.45 0.45]}, ...
    {Lsurf, 'medial',  [0.03 0.03 0.45 0.45]}, ...
    {Rsurf, 'medial',  [0.52 0.03 0.45 0.45]} };

for c = 1:numel(configs)
    cfg = configs{c};
    obj = surface(obj, 'direction', cfg{1}, 'orientation', cfg{2}, 'axes', cfg{3}, rest{:});
end

end

