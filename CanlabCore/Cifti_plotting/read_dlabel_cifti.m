function out = read_dlabel_cifti(fname)

% check for dlabel.nii filetype
if any(regexp(fname,'\.dlabel\.nii$'))

    % check for wb_command on matlab path
    PATH = getenv('PATH');
    if ~contains(PATH,'workbench')
        % add workbench to Matlab Path - only tested on MacOS so far
        if ~ispc
            [~,wbpath] = system('which wb_command');
            if isempty(wbpath) 
                % try default
                setenv('PATH', [PATH ':/Applications/workbench/bin_macosx64']);
            else
                setenv('PATH', [PATH ':' strtrim(wbpath)]);
            end
        end
    end
    
    
    % export dlabel-data to temp textfile
    tmpfile = tempname;
    tmpfile = [tmpfile '.txt'];
    system(['wb_command -cifti-label-export-table ' fname ' 1 ' tmpfile]);
    % read into matlab
    fid  =fopen(tmpfile);
    draw = textscan(fid,'%s\n%d %d %d %d %d');
    fclose(fid);
 
    delete(tmpfile);
    
    % read vertex dat
    out = ft_read_cifti(fname);
    
    % format
    out.label = draw{1};
    out.ID    = single(draw{2});
    out.color = single([draw{3:5}]) / 255;
    out.alpha = single(draw{6}) / 255;
    
else
    error('%s is not a .dlabel.nii-file', fname);
end

end