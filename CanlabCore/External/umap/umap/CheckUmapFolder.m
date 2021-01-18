function ok=CheckUmapFolder(curPath, fileName, doGrandParent)
ok=true;
whereIsFile=fileparts(which(fileName));
if nargin>2 && doGrandParent
    check1=fileparts(curPath);
    [check2, p]=fileparts(whereIsFile);
else
    check1=curPath;
    check2=whereIsFile;
end
if ~isempty(whereIsFile) && ~strcmp(check1, check2)
    if nargin>2 && doGrandParent
        should=fullfile(check1, p);
    else
        should=curPath;
    end
    ok=askYesOrNo(Html.Wrap(['<b>' fileName ...
        '</b> should only be found in '...
        '<br>&nbsp;&nbsp;&nbsp;&nbsp;<b>'...
        should '</b><br><br>...but is also found in'...
        '<br>&nbsp;&nbsp;&nbsp;&nbsp;<b>' whereIsFile ...
        '</b><br><br><i>This could lead to problems, consider first'...
        '<br>clearing your MatLab paths before running run_umap.</i>'...
        '</b><br><br><center><font color="red"><b>'...
        'Continue?</b> </font></center>']), 'Path CONFLICT !!!');
    
end

end
