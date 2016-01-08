function [none] = batch_t_histograms(varargin)
% Creates page(s) of t stat histograms for each subject level contrast in
% set of subject level analyses using image_intensity_histograms.
%
% :Usage:
% ::
%
%     batch_t_histograms([options])
%
% :Optional Inputs:
%
%   **{analysis_dirs}:**
%        run on all contrasts in directories of cell array {analysis_dirs}
%
%       (DEFAULT: use all directories in working directory containing spmT_*.img files)
%
%   'o', 'output_directory':**
%        specify output directory to contain saved .png files

iext = '.img'; %% set up
% image extension


%% parse arguments
dirs = '';
odir = pwd;
i=1;
while i<=numel(varargin)
    if strcmp(varargin{i},'o')
        odir = varargin{i+1};
        i=i+2;
    elseif iscell(varargin{i})
        dirs = [dirs; varargin{i}];
        i=i+1;
    else
        error(['Unrecognized option: ' varargin{i}])
    end
end

if ~exist(odir,'dir'), mkdir(odir); end

%% determine subject level analyis directories and contrast images
% find all contrast files
if isempty(dirs)
    allconfiles = filenames(fullfile('*',['spmT_*' iext]));
else
    allconfiles='';
    for i=1:numel(dirs)
        allconfiles = [allconfiles; filenames(fullfile(dirs{i},['spmT_*' iext]))];
    end
end
% divide by subject level analyses and contrasts
for i=1:numel(allconfiles);
    [subs{i} cons{i}] = fileparts(allconfiles{i});    
end
subs = unique(subs);
cons = unique(cons);


%% make histograms
for c=1:numel(cons)    
    % determine name of output file
    fout = fullfile(odir, [cons{c} '_histograms']);
            
    fprintf('... MAKING histograms for: %s\n',cons{c});    
    % list tmaps
    clear tmaps titles
    n = 1;
    for s=1:numel(subs)
        % determine T map filename
        tmaps{n} = fullfile(subs{s},[cons{c} iext]);
        titles{n} = subs{s};
        n = n+1;
    end
        
    image_intensity_histograms(fout,tmaps,'titles',titles, 'ymax', 15000);
end
