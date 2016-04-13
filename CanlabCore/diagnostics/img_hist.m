function img_hist(imgname,subdir)
% A general function for plotting histograms of any image
% For each subject, comparing across subjects
%
% :Inputs:
%
%   **imgname:**
%        name of image file to make intensity histograms from
%
%   **subdir:**
%        cell array of text strings containing names of individual subject
%        directories (wherein are contained the file specified in imgname 
%        or each subject)
%
% Performs the histogram plot twice, once for CSF space
% and once for gray matter
%
% Locations of gray and CSF masks for each subject must
% be defined in the defaults section of the script.
% (hard-coded)
%
% Start in directory above individual subject results
%
% :Examples:
% ::
%
%    img_hist('beta_0010.img',subdir)
%    img_hist('con_0002.img',{'020827mk' '020829jh' '020903lb'}
%
%    % for batch
%    d = dir('020726ag/beta*img'); d = str2mat(d.name);
%    for i = 1:10:size(d,1)
%        img_hist(deblank(d(i,:)),EXPT.subjects)
%    end
%
%    for i = 2:19, 
%       if i < 10, myz = '000';, else, myz = '00';, end, 
%       img_hist(['con_' myz num2str(i) '.img'],EXPT.subjects);, 
%    end
%
% ..
% Tor Wager
% ..


% ..
%    defaults
% ..
csfname = 'rnnhet1spgr_seg3.img';	% reslice first to space of functionals
grname = 'rnnhet1spgr_seg1.img';
csfpath = '/data/placebo/';		% before ind subject
csfp2 = 'anatomy';							% after ind subject
% dwcard = '02*';									% wildcard defining ind subject directories, e.g., '02*'

% -------------------------------------------------------------------------------------------------

mypwd = pwd;

% get file names

% hP = get_filename(dwcard,['anatomy' filesep imgname]);

for i = 1:length(subdir)

	if i == 1
		hP = fullfile(mypwd,subdir{i},imgname);
		mPcsf = fullfile(csfpath,subdir{i},csfp2,csfname);
		mPg = fullfile(csfpath,subdir{i},csfp2,grname);
	else
		hP1 = fullfile(mypwd,subdir{i},imgname);
		mPcsf1 = fullfile(csfpath,subdir{i},csfp2,csfname);
		mPg1 = fullfile(csfpath,subdir{i},csfp2,grname);

		hP = str2mat(hP,hP1);
		mPcsf = str2mat(mPcsf,mPcsf1);
		mPg = str2mat(mPg,mPg1);

	end
end

% check and remove problematic ones
p1 = hP; p2 = mPcsf;, p3 = mPg;

for i = 1:size(hP,1)

	tmp(i,1) = exist(deblank(hP(i,:)));
	tmp(i,2) = exist(deblank(mPcsf(i,:)));
	tmp(i,3) = exist(deblank(mPg(i,:)));
end

tmp = sum(tmp == 0,2);
tmp = tmp > 0;

hP(tmp,:) = [];
mPcsf(tmp,:) = [];
mPg(tmp,:) = [];

missing = str2mat(subdir{:});
missing = missing(find(tmp),:);

if any(tmp), fprintf(1,'\nMissing required images for\n%s \t ',missing'),end
fprintf(1,'\n')


if isempty(hP), 
	disp('No subjects with all required images. Looking for:'), 
	p1,p2,p3
	return, 
end

csfM = spm_general_hist(hP,mPcsf,'ventricles');
print -P jjon-print1
gM = spm_general_hist(hP,mPg,'graymatter');
print -P jjon-print1

allM = [csfM gM];

diary img_hist.out

disp(['Image ' imgname ' in ' mypwd])
fprintf(1,'\n')
% mean var skew kurt
fprintf(1,'M_csf\tsd_csf\tsk_csf\tkur_csf\tM_gray\tsd_gray\tsk_gray\tkur_gray\t\n');
fprintf(1,'%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t\n',allM');
fprintf(1,'\n')
a = corrcoef(allM);
disp('Correlations')
fprintf(1,'\nM_csf\tsd_csf\tsk_csf\tkur_csf\tM_gray\tsd_gray\tsk_gray\tkur_gray\t\n');
fprintf(1,'%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t\n',a');
fprintf(1,'\n')
diary off


return
