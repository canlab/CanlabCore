function dicom_tarzip(varargin)
% tars all .dcm files in subdirectories you specify, and deletes the
% original .dcm files if successful.
%
% :Inputs:
%
%   none (default): all directories under the current one
%   directory wildcard, e.g., 'run*' to search only those subdirectories.
%   absolute path names should also be ok.
%
% :Examples:
% ::
%
%    dicom_tarzip
%
%    dicom_tarzip('18*')
%
% Or, let's say you have runs within subjects, and dicom files within run dirs.
% Subject dirs all start with NSF*.  run dirs all start with 18*
% Starting from the directory above subjects, you could do this to zip ALL subjects:
% ::
% 
%    subjdirs = dir('NSF*');
%    for s = 1:length(subjdirs)
%        dicom_tarzip([subjdirs(s).name filesep '18*']);
%    end
%
% Needless to say, use this with extreme caution, as it deletes your original files.
% And test it thoroughly with a backup dataset first!!
%
% ..
%    Note: Ash's checksum program would be very useful to add here.
% ..

dirwildcard = '*';

if ~isempty(varargin), dirwildcard = varargin{1}; end

prevwd = pwd;


% list directories

d = dir(dirwildcard);

isd = cat(1, d.isdir);
d = d(isd);
d(strmatch('.', char(d.name))) = [];  % exclude . dirs

if isempty(d)
    fprintf('No matching directories for %s. Quitting.\n', dirwildcard);
    cd(prevwd);
    return
else
    disp('Working on directories:')
    disp(char(d.name))
    disp('  ')
end


basedir = fullfile(pwd, fileparts(dirwildcard));
% % dir returns relative paths, so go there if we have a path name for a
% % wildcard
% if ~isdir(basedir)
%     fprintf('%s is not a directory. Quitting.\n\n', basedir);
%     return
% end
if ~isempty(basedir), cd(basedir); end  % this would be the subject directory if a path name was entered as a wildcard
fprintf('Working in dir: %s\n', basedir)


for i = 1:length(d)
    
    
    
    dicomtarname = [d(i).name filesep d(i).name '_dicoms.tar'];
    
    fprintf('Working on: %s.\n', dicomtarname);
    
    if exist(dicomtarname, 'file')
        fprintf('%s already exists!!! Skipping.\n\n', dicomtarname);
        continue
    end
    
    str = ['!tar cvf ' dicomtarname ' ' d(i).name filesep '*dcm >& tarchecktmp.txt  > tarfilelisttmp.txt'];
    disp(str)
    eval(str)
    errchk = fileread('tarchecktmp.txt');
    isok = isempty(strfind(errchk, 'Error'));
    
    if ~isok
        fprintf('Error in creating: %s .\n', dicomtarname);
        disp(errchk)
        disp('Removing failed tar archive.')
        str = ['!rm ' dicomtarname];
        disp(str)
        eval(str)
        
    else
        % remove original files
        fprintf('%s created successfully. Checking number of files...\n', dicomtarname);
        
        % check length/files
        % TOO LONG str = ['!tar -t ' dicomtarname ' > tarchecktmp.txt'];
        %filelist = fileread('tarfilelisttmp.txt');
        fid = fopen('tarfilelisttmp.txt');
        filelist = textscan(fid, '%s');
        fclose(fid);
        filelist = filelist{1};
        nfiles = length(filelist);
        
        dd = dir([d(i).name filesep '*dcm']);
        nfiles2 = length(dd);
        
        if nfiles == nfiles2
            fprintf('tar archive contains %3.0f files.\n', nfiles);
            fprintf('Removing originals.\n');
            str = ['!rm ' d(i).name filesep '*dcm'];
            disp(str)
            eval(str)
            disp('   ')
        else
            fprintf('tar archive contains %3.0f files, but counted %3.0f that should be there.\n', nfiles, nfiles2);
            disp('Check this. Not removing originals.');
        end
        
    end
    
    !rm tarchecktmp.txt
    !rm tarfilelisttmp.txt
    
    disp('   ')
    
    
    
    
end

if ~isempty(basedir), cd(prevwd); end

end % main function

%dcmf = filenames([d(i).name filesep '*dcm'], 'char');

%%
% dcmf2 = [];
% if ~isempty(dcmf)
%     for ff = 1:size(dcmf, 1)
%        dcmf2 = strcat(dcmf2, [dcmf(ff, :) ' ']);
%     end
%
% end
% str = ['!tar cvf ' d(i).name '_dicoms.tar ' dcmf2];




