function [o2, sig, pcl, ncl] = multi_threshold(dat, varargin)
% Multiple threshold function for statistic_image object for visualization
%
% Usage:
% -------------------------------------------------------------------------
% [o2, sig, pcl, ncl] = multi_threshold(dat, [optional inputs])
%
% Author and copyright information:
% -------------------------------------------------------------------------
%     Copyright (C) 2013  Tor Wager
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
%
% Inputs:
% -------------------------------------------------------------------------
% dat           a statistic_image object
%
% Optional inputs:
% 'poscolors'   followed by cell array of colors for positive values, one per thresh
% 'negcolors'   followed by cell array of colors for negative values, one per thresh
% 'thresh'      followed vector of p-value thresholds, one per thresh
% 'sizethresh'  followed by vector of cluster sizes, one per thresh
%               - this 'prunes' by default, so sizes after first can be 1
%               voxel
%
%               Default thresholds: thresh = [.001 .005 .05], 
%               10 voxels at .001, "pruned"
%
% 'nodisplay'   suppress fmridisplay
%               
% 
%
% Outputs:
% -------------------------------------------------------------------------
% o2            handle to fmridisplay object created by default
% sig           vector of significant voxels at each thresh, for each region
%               - cell array of images in object with matrix of values
%               for each threshold
% pcl           - positive valued clusters cell, one cell per threshold
%               - pass inot mediation_brain_surface_figs.m
% ncl           - positive valued clusters cell, one cell per threshold
%               - pass inot mediation_brain_surface_figs.m
%
% Examples:
% -------------------------------------------------------------------------
%
% [o2, sig, poscl, negcl] = multi_threshold(hr_intercept, 'nodisplay');
% mediation_brain_surface_figs(poscl, negcl);
%
% See also:
% mediation_brain_surface_figs, iimg_multi_threshold,
% mediation_brain_results

% Programmers' notes:
% List dates and changes here, and author of changes



o2 = [];

poscolors = {[1 1 0] [1 .5 0] [.7 0 0]};
negcolors = {[0 0 1] [0 .5 1] [.4 0 .7]};

thresh = [.001 .005 .05];
sizethresh = [10 1 1];

% optional inputs with default values
% -----------------------------------
% - allowable_args is a cell array of argument names
% - avoid spaces, special characters, and names of existing functions
% - variables will be assigned based on these names
%   i.e., if you use an arg named 'cl', a variable called cl will be
%   created in the workspace

allowable_args = {'poscolors', 'negcolors', 'thresh', 'sizethresh', ...
    'nodisplay', 'existingfig'};

default_values = {poscolors, negcolors, thresh, sizethresh, ...
    0, 0};

% define actions for each input
% -----------------------------------
% - cell array with one cell for each allowable argument
% - these have special meanings in the code below
% - allowable actions for inputs in the code below are: 'assign_next_input' or 'flag_on'

actions = {'assign_next_input', 'assign_next_input', 'assign_next_input', 'assign_next_input', ...
    'flag_on', 'flag_on'};

% logical vector and indices of which inputs are text
textargs = cellfun(@ischar, varargin);
whtextargs = find(textargs);

for i = 1:length(allowable_args)
    
    % assign default
    % -------------------------------------------------------------------------
    
    eval([allowable_args{i} ' =  default_values{i};']);
    
    wh = strcmp(allowable_args{i}, varargin(textargs));
    
    if any(wh)
        % Optional argument has been entered
        % -------------------------------------------------------------------------
        
        wh = whtextargs(wh);
        if length(wh) > 1, warning(['input ' allowable_args{i} ' is duplicated.']); end
        
        switch actions{i}
            case 'assign_next_input'
                eval([allowable_args{i} ' = varargin{wh(1) + 1};']);
                
            case 'flag_on'
                eval([allowable_args{i} ' = 1;']);
                
            otherwise
                error(['Coding bug: Illegal action for argument ' allowable_args{i}])
        end
        
    end % argument is input
end

% END DEFAULTS AND INPUTS
% -------------------------------------------------------------------------




sig = {};

nimgs = size(dat.dat, 2); % <- number of 3-D images in dataset

% get significant vectors sig, where sig is cell, one cell per image
for i = 1:length(thresh)
    dat = threshold(dat, thresh(i), 'unc', 'k', sizethresh(i));
    dat = replace_empty(dat);
    for j = 1:nimgs
        sig{j}(:, i) = dat.sig(:, j);
    end
end

% prune
for j = 1:nimgs
    for i = 2:length(thresh)
        sig{j}(:, i) = iimg_cluster_prune(sig{j}(:, i), sig{j}(:, i-1), dat.volInfo);
    end
end


% regions

%r = cell(1, nimgs);
r = cell(1, length(thresh));

for i = 1:length(thresh)
    for j = 1 %:nimgs  % only does first in region anyway
        dat.sig(:, j) = sig{j}(:, i);
        r{i} = region(dat);
    end
    [pcl{i}, ncl{i}] = posneg_separate(reparse_continguous(r{i}));
end

if nodisplay
    return
end


% fmri display
o2 = fmridisplay;
xyz = [-20 -10 -6 -2 0 2 6 10 20]';
xyz(:, 2:3) = 0;

indx = 1;

for j = 1:nimgs
    
    datplot = dat;
    datplot.dat = dat.dat(:, j);
    datplot.p = dat.p(:, j);
    datplot.ste = dat.ste(:, j);
    
    o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 6);
    axh = axes('Position', [0.05 0.4 .1 .5]);
    o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh);
    
    for i = length(thresh):-1:1
        
        datplot.sig = sig{j}(:, i);
        
        mycolors = [negcolors(i) negcolors(i) poscolors(i) poscolors(i)];
        
        o2 = addblobs(o2, region(datplot), 'splitcolor', mycolors, 'wh_montages', indx);
        o2 = addblobs(o2, region(datplot), 'splitcolor', mycolors, 'wh_montages', indx+1);
        
    end
    
    indx = indx + 2;
    
end

end % function
