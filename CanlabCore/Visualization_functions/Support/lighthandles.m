function lh = lighthandles(axishandle)
%lh = lighthandles(axishandle)
% get all light objects associated with the current axis
% see lightRestoreSingle

lh = findobj(axishandle, 'Type', 'Light');

hans = get(axishandle,'Children');
lh2 = find(strcmp(get(hans,'Type'),'light'));
lh = [lh hans(lh2)];

end

