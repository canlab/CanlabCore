function obj = addpoints(obj, xyz, varargin)
% Plots points on fmridisplay objects (e.g., montages of slices)
%
% :Usage:
% ::
%
%     newax = addpoints(obj, xyz, varargin)
%
% Registers handles with the object (referred to as obj)
%
% - enter xyz as n x 3 list of coordinates in mm to plot (world space)
% - Points or text labels or both
% - Flexible slice spacing, colors, marker sizes/styles, axis layout (one row/standard square)
% - axial, saggital, or coronal orientation handled automatically
% - Multiple different sets of points can be plotted in different colors/text labels
%
% :Optional Inputs:
% 
% Takes all inputs of plot_points_on_slice.  See help for additional
% documentation of options.  
%
%   **{'text', 'textcodes'}:**
%        cell array of text values corresponding to points
%
%   **{'condf' 'colorcond'}:**
%        vector of integers to define color conditions
%
%   **'close_enough':**
%        mm within which to plot; defined automatically based on slice distance if not entered
%
%   **'color':**
%        string, 'b', or vector, [1 0 0], to define colors; cell if condf is used, e.g., {'b' 'g'}
%
%   **{'marker', 'MarkerStyle'}:**
%        e.g., 'o', 'v', 's'
%
%   **{'MarkerSize', 'markersize'}:**
%
%   **{'MarkerFaceColor', 'markerfacecolor'}:**
%        see color above
%
% :Examples:
%
% Plot points (i.e., coordinate locations) for xyz coords:
% ::
%
%    o2 = addpoints(o2, DB.xyz, 'MarkerFaceColor', 'b', 'Marker', 'o', 'MarkerSize', 4);
%    o2 = addpoints(o2, DB.xyz, 'text', DB.textcodes, 'condf', DB.condf, 'color', {'b' 'g'});
%    o2 = removepoints(o2);

wh_montage = 1:length(obj.montage); % select which montages; default = all

whm = strcmp(varargin, 'wh_montages') | strcmp(varargin, 'wh_montage') | strcmp(varargin, 'which_montages') | strcmp(varargin, 'which montages');
if any(whm)
    whm = find(whm);
    wh_montage = varargin{whm(1) + 1};
end

% Montages
% -------------------------------------------------------------------------

for i = wh_montage
    
    if ~isfield(obj.montage{i}, 'plotted_point_handles') || isempty(obj.montage{i}.plotted_point_handles)
        
        obj.montage{i}.plotted_point_handles = [];
        
    end
    
    % Plot it
    pointhan = plotpoints(xyz, obj.montage{i}, varargin{:});
    
    obj.montage{i}.plotted_point_handles = [obj.montage{i}.plotted_point_handles pointhan];
    
end



end  % main function





function pointhan = plotpoints(xyz, montagestruct, varargin)


% fixed fields

myview = montagestruct.orientation;
slicemm = montagestruct.slice_mm_coords;
axhan = montagestruct.axis_handles;

% optional inputs

textcodes = [];
texthandles = [];
condf = [];

close_enough = min(diff(slicemm)) ./ 2;
if isempty(close_enough), close_enough = 8; end
    
color = 'k';
marker = 'o';
markersize = 12;
markerfacecolor = 'k';

% ------------------------------------------------------
% parse inputs
% ------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            
            case {'text', 'textcodes'} % do not pass on to slice plot...
                textcodes = varargin{i + 1};
                varargin{i+1} = [];
                varargin{i} = [];
                
            case {'condf' 'colorcond'}, condf = varargin{i + 1};
                
            case 'close_enough', close_enough = varargin{i + 1};
                
            case {'Color','color'} % do not pass on...
                color = varargin{i+1};
                varargin{i+1} = [];
                varargin{i} = [];
                
            case {'marker', 'Marker', 'MarkerStyle'}, marker = varargin{i+1}; varargin{i+1} = [];
                
            case {'MarkerSize', 'markersize'}, markersize = varargin{i+1}; varargin{i+1} = [];
                
            case {'MarkerFaceColor', 'markerfacecolor'}, markerfacecolor = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% SETUP
% -----------------------------------------------

fprintf('Plotting points within %3.2f mm of slices\n', close_enough);

if isempty(condf)
    condf = ones(size(xyz, 1), 1);
    
    if iscell(color)
        error('Color should not be a cell array.');
    end
    color = {color};
    
else
    if ~iscell(color)
        error('When entering condf, enter colors cell with same number of entries.');
    end
end

u = unique(condf);
n = length(u);

% Do the work for each slice
% -----------------------------------------------

pointhan = [];

for i = 1:length(axhan)
    axes(axhan(i));
    
    for j = 1:n % For each color code
        
        % Select coordinates (2D) for plot
        slicexyz = xyz;
        
        switch myview
            case 'axial'
                whcol = 3;
                
            case {'sagg', 'sagittal', 'saggital'}
                
                whcol = 1;
                
            case {'cor', 'coronal'}
                
                whcol = 2;
                
            otherwise
                error('Unknown slice orientation.')
        end
        
        
        
        wh_to_plot = condf == u(j) & abs(slicexyz(:, whcol) - slicemm(i)) <= close_enough;
        
        slicexyz = slicexyz(wh_to_plot, :);
        
        slicexyz(:, whcol) = [];
        
        if ~isempty(slicexyz)
            
            if isempty(textcodes)
                
                pointhan(end + 1) = plot(slicexyz(:, 1), slicexyz(:, 2), '.', 'color', color{j}, 'Marker', marker, 'MarkerSize', markersize, 'MarkerFaceColor', markerfacecolor);
                
            else
                my_text = textcodes(wh_to_plot);
                
                pointhan = [pointhan plottext(slicexyz(:, 1), slicexyz(:, 2), my_text, color{j}, myview)];
                
            end
            
        end
        
        
        
    end % condition color codes
    
end % slices

end % function







function texthandles = plottext(x, y, textcodes, color, orientation)

% adjust for font placement
switch orientation
    case 'sagittal'
        xshift = 6;
    case 'axial'
        xshift = -6;
        
end


texthandles = [];
for j = 1:length(textcodes)
    % text labels
    
    
    if ischar(color)
        
        texthandles(j) = text(x(j) + xshift, y(j), textcodes{j},'Color',color,'FontSize',12,'FontWeight','bold');
        
    else
        texthandles(j) = text(x(j)  + xshift, y(j), textcodes{j},'Color',color,'FontSize',12,'FontWeight','bold');
        
    end
    
end

end % plottext

