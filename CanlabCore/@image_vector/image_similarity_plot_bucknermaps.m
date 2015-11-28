function stats = image_similarity_plot_bucknermaps(obj, varargin)
% Point-biserial correlations between images in fmri_data obj and Bucker
% Lab 7-network maps, with polar plot
%
% :Usage:
% ::
%
%    stats = image_similarity_plot_bucknermaps(obj, 'average');
%
% This is a method for an image_vector object
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2015 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **obj:**
%        An image object with one or more images loaded
%
% :Optional inputs:
%
%   **average:**
%        Calculate average over images in obj with standard errors
%        Useful if obj contains one image per subject and you want
%        to test similarity with maps statistically.
%        Default behavior is to plot each individual image.
%
% :Outputs:
%
%   **stats:**
%        Structure including:
%          - .r, Correlations in [7 networks x images in obj] matrix
%          - .t, T-test (if 'average' is specified)
%          - .line_handles Handles to polar plot lines so you can
%            customize
%          - .fill_handles Handles to polar plot fills so you can
%            customize
%
% :Examples:
% ::
%
%    % corrdat is an fmri_data object with 18 images from searchlight
%    % correlation in it.  Then:
%    stats = image_similarity_plot_bucknermaps(corrdat, 'average');
%
%    % t_diff is a thresholded statistic_image object
%    stats = image_similarity_plot_bucknermaps(t_diff);
%
% :See also:
%
% tor_polar_plot
%
% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

% ..
% DEFAULTS AND INPUTS
% ..

doaverage = 0; % initalize optional variables to default values here.


% optional inputs with default values
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'average', doaverage = 1;
                
                %case 'basistype', basistype = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end



% Load Bucker Lab 1,000FC masks
% ------------------------------------------------------------------------

names = load('Bucknerlab_7clusters_SPMAnat_Other_combined_regionnames.mat');
img = which('rBucknerlab_7clusters_SPMAnat_Other_combined.img');

mask = fmri_data(img);

networknames = names.rnames(1:7);
k = length(networknames);

newmaskdat = zeros(size(mask.dat, 1), k);

for i = 1:k
    
    wh = mask.dat == i;
    
    nvox(1, i) = sum(wh);
    
    newmaskdat(:, i) = double(wh);
    
    
end

mask.dat = newmaskdat;

% Deal with space and empty voxels so they line up
% ------------------------------------------------------------------------

mask = replace_empty(mask); % add zeros back in

mask = resample_space(mask, obj);

obj = replace_empty(obj);

% Correlation
% ------------------------------------------------------------------------

% Point-biserial correlation is same as Pearson's r.
% Gene V. Glass and Kenneth D. Hopkins (1995). Statistical Methods in Education and Psychology (3rd edition ed.). Allyn & Bacon. ISBN 0-205-14212-5.
% Linacre, John (2008). "The Expected Value of a Point-Biserial (or Similar) Correlation". Rasch Measurement Transactions 22 (1): 1154.
% http://www.andrews.edu/~calkins/math/edrm611/edrm13.htm#POINTB

% If both binomial, could use Dice coeff:
%dice_coeff = dice_coeff_image(mask);

% if map or series of maps, point-biserial is better.

% This is done for n images in obj

r = corr(double(obj.dat), double(mask.dat))';

stats.r = r;

if ~doaverage
    
    % Plot values for each image in obj
    [hh, hhfill] = tor_polar_plot({r}, scn_standard_colors(size(r, 2)), {networknames}, 'nonneg');
    
    
elseif doaverage
    
    % Plot mean and se of values
    m = mean(r')';
    se = ste(r')';
    %hh = tor_polar_plot({[max(0, m-se) max(0, m) max(0, m+se)]}, {'r' 'r' 'r'}, {networknames});
    %hh = tor_polar_plot({[m-se m m+se]}, {'r' 'r' 'r'}, {networknames}, 'nonneg');
    
    %hh = tor_polar_plot({m}, {'r'}, {networknames}, 'nonneg');
    [hh, hhfill] = tor_polar_plot({[m+se m m-se]}, {'r' 'r' 'r'}, {networknames}, 'nonneg');
    
    set(hh{1}([1 3]), 'LineWidth', 1); %'LineStyle', ':', 'LineWidth', 2);
    set(hh{1}(2), 'LineWidth', 4);
    set(hhfill{1}([3]), 'FaceAlpha', 1, 'FaceColor', 'w');
    
    %[h, p, ci, stat] = ttest(r');
    [h, p, ci, stat] = ttest(fisherz(r'));
    stats.descrip = 'T-test on Fisher''s r to Z transformed point-biserial correlations';
    stats.networknames = networknames;
    stats.p = p';
    stats.sig = h';
    stats.t = stat.tstat';
    stats.df = stat.df';
        
    disp('Table of correlations');
    disp('--------------------------------------');
    disp(stats.descrip)
    
    print_matrix([m stats.t stats.p stats.sig], {'R_avg' 'T' 'P' 'sig'}, networknames);
    disp(' ');
    
    
end  % doaverage

stats.line_handles = hh;
stats.fill_handles = hhfill;

hhtext = findobj(gcf, 'Type', 'text'); set(hhtext, 'FontSize', 20);
    
end % function

