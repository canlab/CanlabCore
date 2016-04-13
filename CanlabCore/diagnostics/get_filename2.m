function [P,P2,d] = get_filename2(dwcard,varargin)
% :Usage:
% ::
%
%     function [P,P2,d] = get_filename2(search string (as with ls command),[verbose])
%
% Start in directory above individual subject directories
%
% :Inputs:
%
%   **dwcard:**
%        Enter dwcard for wildcard of directories to look in; e.g. '02*'
%
%   **wcard:**
%        Enter wcard for image or file to get - e.g., 'beta_0001.img'
%        This can also include subdirectories (no *'s) 'anatomy/beta*img'
%
% :Outputs:
%
%     Returns list of all files in individual directories  in string matrix
%
%   **P:**
%        file names with full paths
%
%   **P2:**
%        file names only
%
%   **d:**
%        list of directories searched for files
%
% Missing files, or entries in directory that do not contain files, are 
% removed from the list.
% NOT entering a * in dwcard seems to produce an error.
%
% :Examples:
% ::
%
%    P = get_filename2('02*/beta_0001*')
%
% one * is allowed in the directory structure right now, 
% multiple *'s in the filename.
%
% ..
%    Tor Wager 5/2/03
% ..


if length(varargin) > 0, verb = varargin{1};, else, verb = 0;,end

clear d, clear D
 
P = []; P2 = [];
    
[dd,fstr,ext] = fileparts(dwcard);


% get list of directories to search

if any(dd == '*')
    d = dir(dd);
    
    [CWD,newd] = fileparts(dd);
    
    if ~isempty(d), 
        for i = 1:length(d), d(i).name = fullfile(CWD,d(i).name);,end
    end
    
    
    pdname = newd;
    while any(CWD == '*')
        d = [d dir(CWD)];
        
        if ~isempty(d), 
            for i = 1:length(d), d(i).name = fullfile(fileparts(CWD),d(i).name,pdname);,end
        end
        
        [CWD,newd] = fileparts(newd);
        pdname = fullfile(pdname,newd);
        if verb, disp(['Looking in ' CWD]);,end
    end
    
    
    if verb, disp('List of dirs to search:'),str2mat(d.name),disp(' '),end
else
    clear d,
    d(1).name = dd; if isdir(d(1).name), d(1).isdir = 1;, else, d(1).isdir = 0;,end
end

% search for files in dirs

for i = 1:length(d),
	if d(i).isdir, 
		D = fullfile(d(i).name,[fstr ext]);, 
        
        if verb, str = ['dir(''' D ''')'];, disp(str),end
        f = dir(D);
        
        % file name only
        if ~isempty(f), P2 = str2mat(P,str2mat(f(~cat(1,f.isdir)).name));,end
        
        % full path name
        for j = 1:length(f)
            if ~f(j).isdir, 
                if isempty(P)
                    P = fullfile(fileparts(D),f(j).name);
                else
                    P = str2mat(P,fullfile(fileparts(D),f(j).name));,
                end
            end
        end
    end
end

return



