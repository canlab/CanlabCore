function h = qchist(images,Nbins,sparse,XLim,titles)
% This function generates a histogram of activations from a set of
% statistic images.  Generally, you want the images to have a normal
% distribution.  Highly skewed distributions may be indicative of bad
% data.
%
% :Usage:
% ::
%
%     function: h = qchist(dat,Nbins,sparse,XLim)
%
%       This function may generate multiple figuresa with 30 histograms
%       each
%
% :Inputs:
%
%   **images:**
%        List of image file names OR fmri_data object.
%
% :Optional Inputs:
%
%   **Nbins:**
%        Number of bins in each histogram (default = 100)
%
%   **sparse:**
%        flag for generating ONLY histograms (default = 0)
%
%   **XLim:**
%        Xlim (default = [-1 1])
%
%   **titles:**
%        a cell array of subplot titles.  If omitted, titles
%        are inferred from assuming images come from a
%        directory structure that looks like the following:
%        /.../subjname/contrastimage.nii
%
%
%
% ..
%    date: 2/1/2011
%    author: Tor Wager
%
%    Edit: Turned into a function - Scott Schafer 7/14/2012
%    Minor edits - Tor Wager, 8/6/2012
%    Added 'sparse' option - Scott Schafer 1/31/2013
% ..


if nargin < 2
    Nbins = 200;
end

if nargin < 3
    sparse = 0;
end

if ~isa(images, 'image_vector')
    dat = fmri_data(images);
else
    dat = images;
end

dat = remove_empty(dat);
sd = nanstd(dat.dat(:));
m = nanmean(dat.dat(:));

if nargin < 4
    % automatically determine the range of XLim
    % XLim = [-1 1];
    XLim = [m-3.5*sd m+3.5*sd];
end



num_images = size(dat.dat,2);

if num_images > 30, disp('Using separate plots for groups of 30 images.'); end

num_figs = ceil(num_images/30);
h = zeros(num_figs,1);

for j = 1:num_figs
    h = create_figure(['hist' num2str(j)]);
    if num_images >= 30*j
        im_range = 1+30*(j-1):30*j;
    else
        im_range = 1+30*(j-1):num_images;
    end
    for i = im_range
        
        subplot(5,6, i-30*(j-1));
        hold on;
        d = double(dat.dat(:, i));
        d(d == 0 | isnan(d)) = [];
        hist(d, Nbins);
        
        axis tight
        
        try
            hh = findobj(gca, 'Type', 'Patch');
            set(hh, 'EdgeColor', 'none')
        catch
            % maybe old graphics?
        end

        set(gca, 'XLim', XLim);
        
        hh = plot_vertical_line(0);
        set(hh, 'Color', 'r', 'LineWidth', 2);
        
        hh = plot_vertical_line(m-sd);
        set(hh, 'Color', 'r', 'LineStyle', '--');
        
        hh = plot_vertical_line(m+sd);
        set(hh, 'Color', 'r', 'LineStyle', '--');
        
        % figure out subjname, or use passed-in subj names
        if ~exist('titles')
            dd = fileparts(dat.fullpath(i, :));
            [junk, subjname] = fileparts(dd);
        else
            subjname = titles{i};
        end
        title(subjname, 'FontSize', 18)
    end
    
end

ss = get(0, 'ScreenSize');

%Other plots...
if sparse == 0
    plot(dat);
end

end % function

