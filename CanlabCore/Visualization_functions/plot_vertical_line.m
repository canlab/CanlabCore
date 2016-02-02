function han = plot_vertical_line(x,color)
% ..
%    Edited by Tor, Feb 08, to plot multiple lines
% ..

if nargin < 2, color = 'k'; end


ylim = get(gca,'YLim');

for i = 1:length(x)
    han = plot([x(i) x(i)], ylim, color);

end

end
