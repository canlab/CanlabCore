try
    handles.pwd=pwd;cd(handles.fdir);
    [FileName,PathName]=uiputfile('*.mat','Save .mat file...');
    handles.fdir=PathName;cd(handles.pwd);
    varname=inputdlg('Select variable name for cl structure','Select var name',1,{'cl'});
    eval([varname{1} '=handles.cl;'])
    save([PathName FileName],varname{1});
    guidata(hObject,handles);
catch beep,disp(['Error: Operation cancelled or no ''cl'' structure in memory.'])
end