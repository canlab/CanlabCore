function han = makelegend(names, colors, varargin)
% ::
%
%    makelegend(names,colors,[decimal places if numeric entries for names])
%
% :Inputs:
%
%   **names:**
%        must be cell array of names OR a vector of numbers (i.e., thresholds) that will be converted to
%        text
%
%   **colors:**
%        can be cell array of text or rgb values, or matrix of [r g b]
%        values
%
%    Optional Name-Value Pairs:
%     'nofig'     - true to append to existing figure (default: false)
%
%
% :Examples:
% ::
%
%    han = makelegend({'red' 'green' 'blue'}, {'r' 'g' 'b'});
%
%    han = makelegend([.001 .005 .01], {[1 1 0] [1 .5 0] [.7 .3 .3]});


% -------------------------
% Parse optional arguments
% -------------------------
p = inputParser;
addParameter(p, 'nofig', false, @(x) islogical(x) || isnumeric(x));
% addParameter(p, 'decimals', 3, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});
nofig = p.Results.nofig;

% -------------------------
% Convert names to strings if numeric
% Enforce cell array for names, format strings
% -------------------------
if ischar(names)
    names = format_strings_for_legend(names);  % user-defined helper
end

if isnumeric(names)
    % Automatically determine number of decimals
    max_val = max(abs(names));
    if max_val < 0.01
        num_decimals = 4;
    elseif max_val < 0.1
        num_decimals = 3;
    elseif max_val < 10
        num_decimals = 2;
    else
        num_decimals = 0;
    end

    fmt = ['%0.', num2str(num_decimals), 'f'];
    names = arrayfun(@(x) sprintf(fmt, x), names, 'UniformOutput', false);
end

% -------------------------
% Create figure if needed
% -------------------------
if ~nofig
    % create_figure('Legend');
    tor_fig;
end


% if nargin < 3, num_decimals = 3; end

%if nargin < 3, makefig = 1; end
%if makefig, tor_fig;  end

num_entries = length(names);
hold on;

% create_figure('Legend');

% if num_entries < 4
% set(gcf, 'Position', [200 200 200 200]);
% set(gca, 'FontSize', 36);
% else
%     set(gcf, 'Position', [200 200 200 320]);
% set(gca, 'FontSize', 36);
% end

% convert numbers to text
if ~iscell(names) && ~ischar(names)
    
    myfmt = ['%3.' num2str(num_decimals) 'f'];
    
    for i = 1:num_entries
        names2{i} = sprintf(myfmt, names(i));
    end
    
    names = names2;
end
    

for i = 1:num_entries
    
    if iscell(colors)
        mycolor = colors{i};
    else
        mycolor = colors(i,:);
    end
    
    if isstr(mycolor)
        
        % h(i) = plot(0,0, 'ks', 'MarkerFaceColor',mycolor(1), 'MarkerSize', 48, 'Visible', 'off'); %'LineWidth',3);
        h(i) = plot(0,0, 'ks', 'MarkerFaceColor',mycolor(1), 'MarkerSize', 48, 'Visible', 'on'); %'LineWidth',3);
    
    else
        % h(i) = plot(0,0, 'ks', 'MarkerFaceColor',mycolor,  'MarkerSize', 48, 'Visible', 'off'); %'LineWidth',12);
        h(i) = plot(0,0, 'ks', 'MarkerFaceColor',mycolor,  'MarkerSize', 48, 'Visible', 'on'); %'LineWidth',12);


    end
end

han = legend(h,names);



% Set the visibility of plot markers to 'off' after creating the legend
for i = 1:num_entries
    set(h(i), 'Visible', 'off');
end


% if num_entries < 4
% set(han, 'Position', [.25 .35 .45 .45]);
% else
%     set(han, 'Position', [.25 .30 .45 .45]);
% end

axis off

scn_export_papersetup(200);

return
