function Pout = filename_get_new_root_dir(P,newroot,nlevels)
% :Usage:
% ::
%
%     Pout = filename_get_new_root_dir(P,newroot,nlevels)
%
% :Inputs:
%
%   **P:**
%        is files (string matrix)
%
%   **newroot:**
%        Append a new root directory to filenames
%        if newroot is empty, prompt
%
%   **nlevels:**
%        Keeps nlevels dirs deep from old directory:
%          - keeps only filename
%          - keeps filename plus last subdir
%          - keeps filename plus last 2 subdirs, etc.
%
% :Examples:
% ::
%
%    Pout = filename_get_new_root_dir(EXPT.SNPM.P{3},pwd,2)
%
%    % Change root dir to pwd for all images in EXPT
%    for i=1:length(EXPT.SNPM.P), EXPT.SNPM.P{i} = filename_get_new_root_dir(EXPT.SNPM.P{i},pwd,2); end
%
% Note: see regexprep.m for a simpler way!!! This function could be
% improved by using it.
%
% ..
%    tor wager
% ..

if isempty(newroot), newroot = spm_get(-1,'*','Get new root directory'); end

Pout = get_Pout(P(1,:),newroot,nlevels);


for i = 2:size(P,1)
    
    Pout = str2mat(Pout,get_Pout(P(i,:),newroot,nlevels));
    
end

for i = 1:size(Pout, 1)
    isfile(i) = exist(deblank(Pout(i,:)), 'file');
end
fprintf(1,'Found %3.0f valid files.\n', sum(isfile > 0));


return


function Pout = get_Pout(p,str,nlevels);

p = deblank(p);
[d,f,e] = fileparts(p);
filename = [f e];

for i = 2:nlevels
    [d,dd] = fileparts(d);
    filename = fullfile(dd,filename);
end

Pout = fullfile(str,filename);


return

