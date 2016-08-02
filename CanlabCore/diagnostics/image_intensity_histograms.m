function [noreturn] = image_intensity_histograms(fout, imgs, varargin)
%   Makes a sheet (fout.png) of intensity distribution histograms of imgs.
%   Will put an even number of rows of subplots on each page saved.
%   If more than one page is needed, will title outputs fout_1.png, etc.
%
% :Usage:
% ::
%
%     image_intensity_histograms(fout, imgs, [options])
%
% :Inputs:
%
%   **fout:**
%        imfilename to be saved ('.png' will be appended)
%
%   **imgs:**
%        cell array of filenames of images to make histograms of
%
% :Optional Inputs:
%
%   **obj:**
%        im
%
%   **'titles', cellarray:**
%        use strings in cellarray as plot titles (DEFAULT: use file names from imgs)
%
%   **'bins', n:**
%        use n bins in histograms (DEFAULT: 100)
%
%   **'ymax', y:**
%        use y-axis from 0 to y (DEFAULT: 10,000)
%
%   **'xmax', x:**
%        use x-axis from -x to x (DEFAULT: 10)
%
%   **'cols', c:**
%        use c columns of subplots (DEFAULT: 5)
%
%   **'maxrows', r:**
%        use no more than r columns of subplots per page (DEFAULT: 7)
%
%   **'includezeros':**
%        include zero intensities in histograms (DEFAULT: exclude zeros)
%
% ..
% Tor Wager
% ..

%% add path if necessary
if exist('canlab_preproc.m','file') ~= 2
    try
        addpath(genpath('/data/projects/wagerlab/Repository'));
    catch
        error('Failed to add canlab repository /data/projects/wagerlab/Repository to path')
    end
end

%% SET DEFAULTS
% default size of sheet
NCOLS = 5;
NROWS = 7;
% default histogram bins
BINS = 100;
% default xaxis
XMAX = 10;
YMAX = 10000;
% default zero-handling
NOZEROS = true;

%% PARSE INPUTS
% add .png extension to output filename
fout = [regexprep(fout,'\.png$','') '.png'];

% parse options
i=1;
while i<=numel(varargin)
    if ischar(varargin{i})
        switch(varargin{i})
            case {'cols'}
                NCOLS = varargin{i+1};
                i=i+2;
            case {'maxrows'}
                NROWS = varargin{i+1};
                i=i+2;
            case {'bins'}
                BINS = varargin{i+1};
                i=i+2;
            case {'titles'}
                TITLES = varargin{i+1};
                i=i+2;
            case {'xmax'}
                XMAX = varargin{i+1};
                i=i+2;
            case {'ymax'}
                YMAX = varargin{i+1};
                i=i+2;  
            case {'includezeros'}
                NOZEROS = false;
                i=i+1;
            otherwise
                error(['Unrecognized option: ' varargin{i}])
        end
    else
        error(['Unrecognized option: ' varargin{i}])
    end
end
        

%% SETUP
% if titles weren't given, use image filenames
if ~exist('TITLES','var')
    TITLES = imgs;
end

% if need multiple pages, evenly distribute across pages by adjusting NROWS
if numel(imgs) > NCOLS*NROWS
    rows = ceil(numel(imgs)/NCOLS);
    pages = ceil(rows/NROWS);
    NROWS = rows / pages;
end

%% ERROR-CHECKING
if ~exist('fout','var')
    error('no outputfile given')
end

if ~exist('imgs','var')
    error('no images given')
end

if numel(imgs)~=numel(TITLES)
    error('number of titles does not match number of images')
end

%% MAKE PAGE(S) OF HISTOGRAMS
% start figure
figure;
scn_export_papersetup;

%dat=fmri_data(char(imgs));


% add histograms
n = 1;
page = 1;
for i=1:numel(imgs)
    % if sheet is full, start a new one
    if n > NCOLS*NROWS
        fout = regexrep(fout, '.png', ['_' page '.png']);
        fprintf('... SAVING %s\n', fout);
        fprintf('... STARTING NEW SHEET\n');
        saveas(gcf, fout);
        close(gcf);
        n = 1;
        page = page+1;
    end
    
    % determine T map filename
    if exist(imgs{i}, 'file')
        fprintf('... ADDING histogram for %s\n', imgs{i});
        
        % read in T map
        v=spm_vol(imgs{i});
        d=spm_read_vols(v);
%         d = dat.dat(:,i);
        
        % exclude 0s if so desired
        if NOZEROS, d = d(d~=0); end

        % add histogram
        subplot(NROWS, NCOLS, n), hist(d, BINS);
        title(TITLES{i})
        axis([-XMAX XMAX 0 YMAX])
        n = n+1;
    
    else
        warning(['T map not found:' tmap])
    end
end

% save current sheet
if page > 1
    fout = regexprep(fout, '.png',['_' page '.png']);
end
fprintf('... SAVING %s\n\n', fout);
saveas(gcf, fout);
close(gcf);
