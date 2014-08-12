function squeeze_axes(h, squeezepercent)
% squeeze_axes(h, squeezepercent)
%
% Given axis handles h, squeeze them to the left by squeezepercent %.
% This is useful for montages overlapping brain slices.
%
% e.g., for an fmridisplay montage object:
% squeeze_axes(obj.montage{1}.axis_handles, 10)
%
% Negative entries will separate axes, so the operation is reversible:
% squeeze_axes(obj.montage{1}.axis_handles, -10)
%
% you can do this with any axes, e.g.,
% h = findobj(gcf, 'Type', 'axes')
% squeeze_axes(h, 10)

% tor wager, copyright 2011
% 

% (xperc yperc xextent yextent)
% pos = get(h, 'Position');
% pos = cat(1, pos{:});
% 
% spacing = diff(pos(:, 1)); % x-position spacing (usually all the same)
% 
% shiftby = (spacing .* squeezepercent ./ 100); 
% 
% shiftby = [0; shiftby .* (1:length(h)-1)']; % cumulative
% 
% pos(:, 1) = pos(:, 1) - shiftby;
% 
% for i = 1:length(h)
% set(h(i), 'Position', pos(i, :));
% end
% 
% % left bottom justify
pos = get(h, 'Position');
pos = cat(1, pos{:});
pos(:, 1) = pos(:, 1) - min(pos(:, 1)) + .05;
pos(:, 2) = pos(:, 2) - min(pos(:, 2)) + .05;

for i = 1:length(h)
set(h(i), 'Position', pos(i, :));
end

% % shrink figure
% % ----------------------------------------------------
% % [lastpos, wh] = max(pos(:, 1));
% % lastpos = lastpos + pos(wh, 2) + .05; % right-most value, in %
% 
% f = get(h(1), 'Parent');
% fpos = get(f, 'Position');
% % (xstart ystart from bottom L, xstretch ystretch)
% 
% %fpos(3) = fpos(3) .* lastpos;
% fpos(3) = fpos(3) - (fpos(3) .* squeezepercent/100);
% 
% set(f, 'Position', fpos);
% 
% % now make x and y axis positions larger
% % ----------------------------------------------------
% pos = get(h, 'Position');
% pos = cat(1, pos{:});
% 
% pos(:, [3 4]) = pos(:, [3 4]) + (pos(:, [3 4]) .* squeezepercent/100);
% 
% for i = 1:length(h)
% set(h(i), 'Position', pos(i, :));
% end

end % function