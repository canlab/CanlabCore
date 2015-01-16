function lh = lighthandles(axishandle)
%lh = lighthandles(axishandle)
% get all light objects associated with the current axis
% see lightRestoreSingle

hans = get(axishandle,'Children');
lh = find(strcmp(get(hans,'Type'),'light'));
lh = hans(lh);

return

