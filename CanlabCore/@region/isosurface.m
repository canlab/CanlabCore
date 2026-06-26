function [surface_handles, colors, patch_cell] = isosurface(r, varargin)
% isosurface Create a series of surfaces in different colors, one per region.
%
% Render each element of a region object array as a 3-D isosurface using
% imageCluster, with options for a single color across regions and for
% matching colors across left/right hemispheres.
%
% :Usage:
% ::
%
%     [surface_handles, colors, patch_cell] = isosurface(r, [optional arguments])
%
% :Inputs:
%
%   **r:**
%        A region-class object array.
%
% :Optional Inputs:
%
%   **Any optional inputs to imageCluster:**
%        e.g., 'alpha' followed by transparency value, 'sd', or
%        'smoothbox'.
%
%   **'colors':**
%        Followed by a single color in { } or a cell array of multiple
%        colors.
%
%   **'nomatchleftright' or 'nosymmetric':**
%        Do not match colors across hemispheres (left/right). The
%        default matches L/R and may override user-supplied colors.
%
% :Outputs:
%
%   **surface_handles:**
%        Vector of patch handles, one per region (legacy form; coerces
%        Patch objects to doubles).
%
%   **colors:**
%        Cell array of [r g b] color triplets used for each region.
%
%   **patch_cell:**
%        Cell array of patch object handles, one per region.
%
% :Examples:
% ::
%
%     atlasfile = which('Morel_thalamus_atlas_object.mat');
%     load(atlasfile)
%
%     surface_handles = isosurface(r);
%     surface_handles = isosurface(r, 'alpha', .5);
%     surface_handles = isosurface(r, 'alpha', .5, 'nomatchleftright');
%
%     view(135, 30);
%     lightRestoreSingle;
%     lightFollowView;
%
%     load(which('CIT168_MNI_subcortical_atlas_object.mat'));
%     r = atlas2region(atlas_obj);
%
%     surface_handles = isosurface(r, 'alpha', .5, 'color', {[.3 .6 .4] [.5 .4 .2]});
%     surface_handles = isosurface(r, 'alpha', .5, 'color', {[.3 .6 .4] [.5 .4 .2]}, 'nomatchleftright');
%     p = addbrain('hires right');
%     lightFollowView;
%
% :See also:
%   - atlas/isosurface
%   - imageCluster
%   - match_colors_left_right


k = length(r);

% ..
%    DEFAULTS AND INPUTS
% ..

% input_args = varargin;           % save for later
colors = scn_standard_colors(k); % generate colors
matchcolorsleftright = true;

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            % do nothing for inputs passed on to imageCluster
            case {'alpha'}  
                
            case {'color', 'colors'}
                colors = varargin{i+1}; varargin{i+1} = []; varargin{i} = [];

                if ~iscell(colors), colors = {colors}; end
                
                if length(colors) == 1 % single color; matchcolorsleftright will overwrite this
                    matchcolorsleftright = false;
                end
                
            case {'nomatchleftright', 'nosymmetric'}, matchcolorsleftright = false; varargin{i} = [];
                
            case {'sd' 'smoothbox'}
                % do nothing, but pass these into imageCluster later
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if iscolumn(colors), colors = colors'; end

% Handle single-color input, and other cases where too few colors:
while length(colors) < k, colors = [colors colors]; end

if matchcolorsleftright
    colors = match_colors_left_right(r);
end

colors = colors(1:length(r));

% -------------------------------------------------------------------------
% Make surfaces
% -------------------------------------------------------------------------

cl = region2struct(r);

surface_handles = []; % This appears to be a bug, as it coerces the Patch object into a double. 05/23/2023 MS
patch_cell = {};

for i = 1:k
    
    try
        out = imageCluster('cluster', cl(i), 'color', colors{i}, varargin{:});
        
        surface_handles(i) = out; % This appears to be a bug, as it
%         coerces the Patch object into a double. 05/23/2023 MS
        patch_cell{i} = out;
    catch
        disp('Error imaging isosurface; too few voxels?')
        surface_handles(i) = [];
        patch_cell{i}=[];
    end
    
%     % Isocaps, if needed
%     
%     % pp = [];
%     isocap = isocaps(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, V, mythresh);
%     p(end + 1) = patch(isocap, 'FaceColor', 'interp','EdgeColor', 'none', 'FaceAlpha',1);


end


% -------------------------------------------------------------------------
% Lighting, etc.
% -------------------------------------------------------------------------

lightRestoreSingle;
lighting gouraud;
axis vis3d image
material dull

end % function


