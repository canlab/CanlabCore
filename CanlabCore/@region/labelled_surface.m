function [centroid, p] = labelled_surface(r, varargin)
% labelled_surface Plot labelled surfaces with centroids and text annotations.
%
% Render a region object as transparent isosurfaces on top of a brain
% surface, label each region with its short title and a sequential
% number, and return per-region centroid coordinates and the underlying
% brain surface handle. The text label is offset toward the camera if
% 'popout' is requested.
%
% :Usage:
% ::
%
%     [centroid, p] = labelled_surface(r, [optional inputs])
%
% :Inputs:
%
%   **r:**
%        A region-class object array whose elements will be rendered.
%
% :Optional Inputs:
%
%   **'surface_keyword':**
%        Keyword passed to addbrain to draw the underlying brain
%        surface. Default: 'transparent_surface'.
%
%   **'colormap':**
%        Function handle returning an [n x 3] colormap. Default:
%        @colorcube.
%
%   **'font_arguments':**
%        Cell array of additional name/value pairs forwarded to text(),
%        e.g., {'FontSize', 12, 'FontWeight', 'bold'}.
%
%   **'popout':**
%        Logical; if true, text labels are pushed toward the camera so
%        they are visible above the surface. Default: false.
%
% :Outputs:
%
%   **centroid:**
%        Cell array of centroid coordinates for each rendered surface.
%
%   **p:**
%        Handle to the underlying brain surface plot.
%
% :Examples:
% ::
%
%     [centroid, p] = labelled_surface(r, 'surface_keyword', 'left_cutaway', ...
%         'colormap', @jet, 'font_arguments', ...
%         {'FontSize', 12, 'FontWeight', 'bold'}, 'popout', true);
%
% :See also:
%   - addbrain
%   - isosurface
%   - format_strings_for_legend
%
% ..
%    Author: Michael Sun, Ph.D. 05/16/2024
% ..

% Parse optional inputs
parser = inputParser;
addParameter(parser, 'surface_keyword', 'transparent_surface');
addParameter(parser, 'colormap', @colorcube);
addParameter(parser, 'font_arguments', {});
addParameter(parser, 'popout', false);
parse(parser, varargin{:});
surface_keyword = parser.Results.surface_keyword;
colormap_func = parser.Results.colormap;
user_font_arguments = parser.Results.font_arguments;
popout = parser.Results.popout;


% Create figure
% figure;

% Add brain surface
p = addbrain(surface_keyword);
set(p, 'FaceAlpha', 0.2);
lightFollowView;

% Initialize centroid cell array
centroid = cell(numel(r), 1);

% Add text and plot surfaces
counter = 1;
cmap = colormap_func(numel(r));
for i = 1:numel(r)
    try
        [~,~,surface_handles] = isosurface(r(i), 'nomatchleftright', 'color', cmap(i,:));
        centroid{i} = mean(surface_handles{1}.Vertices);
    
        % Format string for legend
        if ~isempty(r)
            label_text = format_strings_for_legend([num2str(counter), '. ', r(i).shorttitle]);
        else
            label_text = format_strings_for_legend([num2str(counter)]);
        end
    
        % Default font arguments including color mapping
        default_font_arguments = {'FontSize', 12, 'FontWeight', 'bold', 'Color', cmap(i,:)};
    
        % Combine default and user-provided font arguments
        font_arguments = combineFontArguments(default_font_arguments, user_font_arguments);
        
        set(surface_handles{1}, 'FaceAlpha', 0.5);
        % set(surface_handles{1}, 'FaceAlpha', 0.0);
    
        % Adjust the text position if popout is enabled
        % Very hacky way to do this.
        text_coords = centroid{i};
        if popout
            % Get current view direction vector
            cam_pos = get(gca, 'CameraPosition');
            cam_tgt = get(gca, 'CameraTarget');
            view_dir = cam_tgt - cam_pos;
            view_dir = view_dir / norm(view_dir); % Normalize the direction vector
            
            % Use a fixed significant offset to move text to the front in the view direction
            offset = -1000; % Adjust this value as needed to ensure text visibility
            text_coords = text_coords + view_dir * offset;
        end
        
        
        text(text_coords(1), text_coords(2), text_coords(3), label_text, font_arguments{:});
    catch

    end
    counter = counter + 1;
end

% Draw and snapshot
drawnow;
snapnow;

end

function combined_args = combineFontArguments(defaults, user_args)
    % Combine default and user-provided font arguments, with user arguments taking precedence
    default_struct = struct(defaults{:});
    user_struct = struct(user_args{:});
    
    fields = fieldnames(user_struct);
    for i = 1:numel(fields)
        default_struct.(fields{i}) = user_struct.(fields{i});
    end
    
    combined_args = reshape([fieldnames(default_struct) struct2cell(default_struct)]', 1, []);
end
