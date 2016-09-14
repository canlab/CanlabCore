function han = makelegend(names, colors, makefig)
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
% :Examples:
% ::
%
%    han = makelegend({'red' 'green' 'blue'}, {'r' 'g' 'b'});
%
%    han = makelegend([.001 .005 .01], {[1 1 0] [1 .5 0] [.7 .3 .3]});

% Enforce cell array for names, format strings
if ischar(names)
    names = format_strings_for_legend(names);
end

if nargin < 3, num_decimals = 3; end

%if nargin < 3, makefig = 1; end
%if makefig, tor_fig;  end

num_entries = length(names);

create_figure('Legend');

if num_entries < 4
set(gcf, 'Position', [200 200 200 200]);
set(gca, 'FontSize', 36);
else
    set(gcf, 'Position', [200 200 200 320]);
set(gca, 'FontSize', 36);
end

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
        
        h(i) = plot(0,0, 'ks', 'MarkerFaceColor',mycolor(1), 'MarkerSize', 48); %'LineWidth',3);
    
    else
        h(i) = plot(0,0, 'ks', 'MarkerFaceColor',mycolor,  'MarkerSize', 48); %'LineWidth',12);
        
    end
end

han = legend(h,names);

if num_entries < 4
set(han, 'Position', [.25 .35 .45 .45]);
else
    set(han, 'Position', [.25 .30 .45 .45]);
end

axis off

scn_export_papersetup(200);

return
