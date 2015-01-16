function EXPT = getfunctnames(EXPT,varargin)
% EXPT = getfunctnames(EXPT,[image wildcard],[save in field],[search in subjects],[startdir])
%
% Examples:
% 
% for contrast/FIR model images 
% EXPT = getfunctnames2(EXPT,'dx*img','FIR.dxbetas');
% EXPT = getfunctnames2(EXPT,'swa*img','FILES.im_files_neg','functional/neg_scans/scan*');
%
% Example to get a set of contrast images and store them in EXPT.SNPM.P
% for i = 1:9
% EXPT.tmp = {};
% EXPT = getfunctnames2(EXPT,['con_000' num2str(i) '.img'],'tmp'); EXPT.tmp = str2mat(EXPT.tmp{:});
% EXPT.SNPM.P{i} = EXPT.tmp;
% end

wcard = 'w*img';   if length(varargin) > 0, wcard = varargin{1}; end
fldname = 'nravols';   if length(varargin) > 1, fldname = varargin{2}; end
subjects = '.';   if length(varargin) > 2, subjects = varargin{3}; end
startdir = pwd;   if length(varargin) > 3, startdir = varargin{4}; end

origdir = pwd;
if isdir(startdir), cd(startdir), else, error('Start directory is invalid--not a directory.'); end

fprintf(1,'Finding %s images in subjects %s, storing in %s\n',wcard,subjects,fldname); 

if ~isfield(EXPT, 'subjects')
    disp('You must have a field in EXPT called EXPT.subjects that lists all subject names.')
    disp('Each subject name should go in a cell in EXPT.subjects.');
    error('Exiting now.');
end

% simpler, faster way
for s = 1:length(EXPT.subjects)

    fprintf(1,'%s.',EXPT.subjects{s})

    dd = [];   	p = pwd;

    %for i = 1:4 
        
        %d = dir(fullfile(EXPT.subjects{s} '/scan' num2str(i) '/' wcard]); d = str2mat(d.name); 
	
        %d = get_filename2(fullfile(EXPT.subjects{s},subjects,wcard));
        d = filenames(fullfile(EXPT.subjects{s},subjects,wcard),'char');
        
        %dp = fullfile(p,EXPT.subjects{s},['scan' num2str(i)],filesep);
        
        % add current path to start of names for full paths
        dp = fullfile(p,filesep);
        dp = repmat(dp,size(d,1),1);
        d  = [dp d];

        dd = [dd; d];
%end

str = (['EXPT.' fldname '{s} = dd;']);
eval(str)
end

d2 = d(:,end-7:end-4);
if length(unique(d2,'rows')) < length(d2)
    disp('Multiple files with same number!  Sorting off -- check output lists for order correctness.')
    
else
        

    fprintf(1,'\nChk sorting. ');

    str = (['EXPT.' fldname ' = sort_image_filenames([EXPT.' fldname ']);']);
    disp(str)
    eval(str)
end

cd(origdir)

return

