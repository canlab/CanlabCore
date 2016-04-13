function P = get_filename(dwcard,wcard,varargin)
% :Usage:
% ::
%
%     function P = get_filename(dir_wcard,file_wcard,[verbose])
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
% :Output:
%
%     Returns list of all files in individual directories  in string matrix
%
% Missing files, or entries in directory that do not contain files, are 
% removed from the list.
%
% NOT entering a * in dwcard seems to produce an error.
%
% :Examples:
% ::
%
%    P = get_filename('02*','beta_0001*')
%    P = get_filename('02*','beta_0001.img')
%    P = get_filename('02*','anatomy/nnhe*_seg1.img')
%    P = get_filename('020515sp*','anatomy/nnhe*_seg1.img')
%
% one * is allowed in first string, multiple *'s in second,
% as long as they are in the filename, not directory names!
%
% ..
%    Tor Wager 10/8/02
% ..

if length(varargin) > 0, verb = varargin{1};, else, verb = 0;,end

clear d, clear D
d = dir(dwcard); dind = 1;, 
for i = 1:length(d),
	if d(i).isdir, 
		D{dind} = d(i).name;, 
		dind = dind + 1;, 
	end, 
end

[ld f e] = fileparts(wcard);
if isempty(ld)
	P = tor_list_files(D,wcard);
	
else
	for i = 1:length(D)
		D{i} = [D{i} filesep ld];
	end
	P = tor_list_files(D,[f e]);
end
P = str2mat(P);

for i = 1:size(P,1), 
	if isempty(deblank(P(i,:))), 
		%warning(['No file found for ' D{i}]), 
		pp(i) = 1;, 
	else, 
		pp(i) = 0;, 
	end
end


for i = 1:size(P,1)
    if isempty(deblank(P(i,:))), pp(i) = 1;, end
end

P = P(~pp,:);

return



