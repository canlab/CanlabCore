function [centroid, p] = labelled_surface(r, varargin)
% labelled_surface - Plots labelled surfaces with centroids and text annotations
% 
% Syntax:  [centroid, p] = labelled_surface(r, varargin)
%
% Inputs:
%    r - Cell array of surfaces to be plotted
%    varargin - Additional parameters for customization (e.g., surface_keyword, colormap, font_arguments)
%
% Outputs:
%    centroid - Cell array of centroid coordinates for each surface
%    p - Handle to the brain surface plot
%
% Example: 
%    [centroid, p] = labelled_surface(r, 'surface_keyword', 'left_cutaway', 'colormap', @jet, 'font_arguments', {'FontSize', 12, 'FontWeight', 'bold'}, 'popout', true);
%
% See also: addbrain, isosurface, format_strings_for_legend

% Author: Michael Sun, Ph.D. 05/16/2024

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
figure;

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
