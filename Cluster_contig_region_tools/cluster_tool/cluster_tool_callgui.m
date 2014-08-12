function out=cluster_tool_callgui(guitocall)

% returns handles.data from guitocall and closes the gui. note that
% guitocall must use uiwait and uiresume WITHOUT the need to close.

eval(['h=' guitocall ';'])
d=guidata(h);
out=d.data;
close(h);