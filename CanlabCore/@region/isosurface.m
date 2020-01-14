function surface_handles = isosurface(r, varargin)
% Create a series of surfaces in different colors, one for each region
% - Options for single color
%
% surface_handles = isosurface(r, [optional arguments])
%
% optional arguments: 
% - Any optional inputs to imageCluster, e.g., 'alpha'
% - 'colors', followed by single color in { } or cell array of multiple colors
% - 'nomatchleftright' or 'nosymmetric', do not match colors across hemispheres (left/right)
%               Note: The default matches, and may override your colors.
%
% Examples:
% atlasfile = which('Morel_thalamus_atlas_object.mat');
% load(atlasfile)
%
% surface_handles = isosurface(r);
% surface_handles = isosurface(r, 'alpha', .5);
% surface_handles = isosurface(r, 'alpha', .5, 'nomatchleftright');
%
% view(135, 30);
% lightRestoreSingle;
% lightFollowView;
%
% load(which('CIT168_MNI_subcortical_atlas_object.mat'));
% r = atlas2region(atlas_obj);
%
% surface_handles = isosurface(r, 'alpha', .5, 'color', {[.3 .6 .4] [.5 .4 .2]});
% surface_handles = isosurface(r, 'alpha', .5, 'color', {[.3 .6 .4] [.5 .4 .2]}, 'nomatchleftright');
% p = addbrain('hires right');
% lightFollowView;


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
                if length(colors) == 1 % single color; matchcolorsleftright will overwrite this
                    matchcolorsleftright = false;
                end
                
            case {'nomatchleftright', 'nosymmetric'}, matchcolorsleftright = false; varargin{i} = [];
                
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


% -------------------------------------------------------------------------
% Make surfaces
% -------------------------------------------------------------------------

cl = region2struct(r);

surface_handles = [];

for i = 1:k
    
    out = imageCluster('cluster', cl(i), 'color', colors{i}, varargin{:});
    
    surface_handles(i) = out;
    
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


