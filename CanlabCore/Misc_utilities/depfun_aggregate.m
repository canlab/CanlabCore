function dependencies = depfun_aggregate(funname, excludes, save_dir)
% depfun_aggregate is for packaging all the files needed for a function in one spot
%
% :Usage:
% ::
%
%     dependencies = depfun_aggregate(funname, excludes, save_dir)
%
% :Inputs:
%
%   **funname:**
%        function to analyze
%
%   **excludes:**
%        cellstrs to exclude from the matches
%
%   **save_dir:**
%        directory to copy files to - defaults to 'tmp_files'
%
% :Example:
% ::
%
%    dependencies = depfun_aggregate('whole_brain_fir', {'spm2', 'spm5'}, 'fir_files')


    dependencies = depfun(funname, '-quiet');
    dependencies = remove_str_from_cell(dependencies, '/MATLAB');
    
    if(~isempty(excludes))
        excludes = cellstr(excludes);
        for i=1:length(excludes)
            dependencies = remove_str_from_cell(dependencies, excludes{i});
        end
    end
    
    if(~exist('save_dir', 'var') || isempty(save_dir))
        save_dir = 'tmp_files';
    end
    mkdir(save_dir);
    fprintf(1, 'Saving to %s\n', save_dir);

    %NB: copyfile() should be used for platform-independence, but 
    %its handling permissions incorrectly on my machine - Matthew!rm 
    for i = 1:length(dependencies)
        system(['cp ' dependencies{i} ' ' save_dir '/']);
%         copyfile(dependencies{i}, 'tmp_files/'); 
        fprintf(1, 'Copied %s to %s.\n', dependencies{i}, save_dir);
    end
end
