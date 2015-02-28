function plot(d)

x = get_x(d);
y = get_y(d);
labels=[];
if(size(y,2)==1)
    uy = unique(y);
    nc = length(uy);
    labels=1:nc;
else
    nc = size(y,2);
    uy = 1:nc;
end
    

cm = jet(max(length(uy), 32));
cml = size(cm, 1);
markers = {'square', 'o', 'x', '+', '^', '*', 'diamond', '<', '>', 'pentagram', 'v', 'hexagram'};
mi = 0;
if nc > length(markers), markers = {'o'}; end

washeld = ishold;
for i = 1:nc
	cmi = round(1+(cml-1)*(i-1)/(nc-1));
	colour = cm(cmi, :);
	mi = 1 + rem(mi, length(markers));
	marker = markers{mi};
	h(i) = plot(x(y==uy(i)), 1), x(y==uy(i), 2), 'color', colour, ...
                'marker', marker, 'markerfacecolor', colour, 'markeredgecolor', colour, 
            'linestyle', 'none','markersize',6);
	hold on
end
if ~washeld, hold off, end
figure(gcf)
drawnow

