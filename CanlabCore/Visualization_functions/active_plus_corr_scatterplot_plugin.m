disp('Running: active_plus_corr_scatterplot_plugin');

% ------------------------------------------------
% define some defaults
% ------------------------------------------------

%if length(varargin) > 0, overlay = varargin{1};, end
if ~isfield(EXPT,'overlay') || isempty(EXPT.overlay),
    EXPT.overlay = which('scalped_single_subj_T1.img');
end
if isempty(which(EXPT.overlay)), warning('Cannot find overlay image on path!'); end
disp(['Overlay: ' EXPT.overlay]);

clusterfield = 'timeseries';

if ~(exist('EXPT') == 1) 
    error('You must load a valid EXPT structure with required fields to run.  Best to run this function from the GUI.');
end

% ------------------------------------------------
% get a directory and go to it.
% ------------------------------------------------
dirname = spm_get(-1,'rob*','Choose robust* or robseed* directory.');
if dirname(end) == '.', dirname = dirname(1:end-2); elseif dirname(end) == filesep, dirname = dirname(1:end-1); end

disp(['Going to: ' dirname])
cd(dirname)

% ------------------------------------------------
% tell whether we're in a robust* or robseed* directory, and do the
% appropriate thing.
% ------------------------------------------------
[dd,robdir] = fileparts(dirname);

dirtype = robdir(1:end-4);          % robust or robseed
% dirnum = str2num(robdir(end-3:end)); % which number -- to get correct data
dirnum = find(EXPT.SNPM.connums == str2num(robdir(end-3:end))); % which number -- to get correct data

% ------------------------------------------------
% get clusters from results image
% save average data in 'timeseries' field, cl.timeseries
% ------------------------------------------------
cl = image2clusters(EXPT.overlay);              % could be masked ordered or not.
cl = tor_extract_rois(EXPT.SNPM.P{dirnum},cl);  

disp('Based on directory name, extracted cluster average data from: '); disp(EXPT.SNPM.P{dirnum});
disp('Stored data from images in cl.timeseries');

switch dirtype
    
    case 'robust'

        disp(EXPT)
        disp('Note: behavioral/other covariates should be entered in EXPT.cov.');
        whcol = spm_input('Enter column of EXPT.cov to plot against data: ',[],[],1);
        
        for i = 1:length(cl), cl(i).xdat = EXPT.cov(:,whcol); end 
        
        docovs = 0;
        if isfield(EXPT,'cov') && size(EXPT.cov,2) > 1
            docovs = spm_input('Adjust for other columns of EXPT.cov? (1/0) ',[],[],1);
        end
        if docovs, for i = 1:length(cl), cl(i).nuisance = EXPT.cov;  cl(i).nuisance(:,whcol) = []; end , end
        
        
    case 'robseed'
        
        disp(EXPT)
        disp('Note: seed data should be entered in EXPT.seeds.');
        if isfield(EXPT,'seednames'), disp(str2mat(EXPT.seednames)); end
        whcol = spm_input('Enter column of EXPT.seeds to plot against data: ',[],[],1);
        
        for i = 1:length(cl), cl(i).xdat = EXPT.seeds(:,whcol); end 
        
        
    otherwise
        error('You must choose a robust or robseed directory.');
        
end



% now we're ready to make the plot, and set the button up function to the
% command below:
[xdat,ydat] = cluster_interactive_scatterplot(cl);

