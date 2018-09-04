function [o2, dat, sig, pcl, ncl] = multi_threshold(dat, varargin)
% Multiple threshold function for statistic_image object for visualization
%
% :Usage:
% ::
%
%    [o2, dat, sig, pcl, ncl] = multi_threshold(dat, [optional inputs])
%
% ..
%     Author and copyright information:
%
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
% ..
%
% :Inputs:
%
%   **dat:**
%        a statistic_image object
%
% :Optional Inputs:
%
%   **poscolors:**
%        followed by cell array of colors for positive values, one per thresh
%
%   **negcolors:**
%        followed by cell array of colors for negative values, one per thresh
%
%   **thresh:**
%        followed vector of p-value thresholds, one per thresh
%
%   **sizethresh:**
%        followed by vector of cluster sizes, one per thresh
%           - this 'prunes' by default, so sizes after first can be 1
%             voxel
%           - Default thresholds: thresh = [.001 .005 .05],
%             10 voxels at .001, "pruned"
%
%   **nodisplay:**
%        suppress fmridisplay
%
%   **o2:**
%       followed by an existing fmridisplay object
%           - will remove blobs and re-use montages
%       Or enter any fmridisplay object, and it will be detected.
%
%   **writestats**
%        write out thresholded images
%
%   **noplot**
%        don't plot anything (for writing stats out)
%
%   **wh_montages** Which existing montages in o2 fmridisplay object to
%   plot on. Optional. Keyword followed by vector of integers for which.
%
% :Outputs:
%
%   **o2:**
%        handle to fmridisplay object created by default
%
%   **sig:**
%        vector of significant voxels at each thresh, for each region
%               - cell array of images in object with matrix of values
%               for each threshold
%
%   **dat:** 
%        thresholded statistic_image object
%
%   **pcl:**
%        positive valued clusters cell, one cell per threshold
%           - FIRST image in object only
%           - pass into mediation_brain_surface_figs.m
%
%   **ncl:**
%        positive valued clusters cell, one cell per threshold
%           - FIRST image in object only
%           - pass into mediation_brain_surface_figs.m
%
% :Examples:
% ::
% % Example 1: A complete analysis of a sample dataset
% % ----------------------------------------------------------------------
% % Load a sample dataset and do a group t-test
% img_obj = load_image_set('emotionreg');         % Load a dataset
% t = ttest(img_obj, .005, 'unc');                % Do a group t-test
% 
% [o2, sig, poscl, negcl] = multi_threshold(t);
% 
% % Make a table of the regions at the highest threshold
% r = table(poscl{1});                            % r is a region object, now with labels for regions attached
% 
% % Show each individual region significant at the highest threshold
% montage(r, 'colormap', 'regioncenters');
%
% % Multi-threshold with a custom threshold.  Do not display immediately,
% % but use mediation_brain_surface_figs to display results on surfaces instead:
%
%    [o2, sig, poscl, negcl] = multi_threshold(t, 'nodisplay', 'thresh', [.001 .005 .01], 'k', [10 1 1]);
%    mediation_brain_surface_figs(poscl, negcl);
%
%    % Create empty montage set and (re)use it:
%    o2 = canlab_results_fmridisplay([], 'compact2', 'noverbose');
%    o2 = multi_threshold(out.t, 'o2', o2);
%
% :See also:
% mediation_brain_surface_figs, iimg_multi_threshold,
% mediation_brain_results, statistic_image.threshold


% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%    Tor Wager - 2018 July - update to take advantage of new object display
%    code, and add documentation
%    Tor - Aug 2018 - can now pass in wh_montages and other optional args
%    to use existing montages. Pass out dat as well.
% ..

poscolors = {[1 1 0] [1 .5 0] [.7 0 0]};
negcolors = {[0 0 1] [0 .5 1] [.4 0 .7]};

%o2 = [];

thresh = [.001 .005 .05];
sizethresh = [10 1 1];

c=0;str_args={};
for a=1:length(varargin)
    if isstr(varargin{a})
        c=c+1;
        str_args{c}=varargin{a};
    end
end

doWrite=0;
doPlot=1;

for s=1:length(str_args)
    
    if ~isempty(strfind(str_args{s},'writestats'))
        doWrite=doWrite+1;
    else
        %         doWrite=0;
    end
    
    
    if ~isempty(strfind(str_args{s},'noplot'))
        doPlot=doPlot-1;
    else
        %         doPlot=1;
    end
end


% optional inputs with default values
% -----------------------------------
% - allowable_args is a cell array of argument names
% - avoid spaces, special characters, and names of existing functions
% - variables will be assigned based on these names
%   i.e., if you use an arg named 'cl', a variable called cl will be
%   created in the workspace

allowable_args = {'poscolors', 'negcolors', 'thresh', 'sizethresh', ...
    'nodisplay', 'existingfig', 'o2'};

default_values = {poscolors, negcolors, thresh, sizethresh, ...
    0, 0, []};

% define actions for each input
% -----------------------------------
% - cell array with one cell for each allowable argument
% - these have special meanings in the code below
% - allowable actions for inputs in the code below are: 'assign_next_input' or 'flag_on'

actions = {'assign_next_input', 'assign_next_input', 'assign_next_input', 'assign_next_input', ...
    'flag_on', 'flag_on', 'assign_next_input'};

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

% special for o2 : detect

wh = cellfun(@(x) isa(x, 'fmridisplay'), varargin);
if any(wh), o2 = varargin{wh}; end


% END DEFAULTS AND INPUTS
% -------------------------------------------------------------------------




sig = {};

nimgs = size(dat.dat, 2); % <- number of 3-D images in dataset

% get significant vectors sig, where sig is cell, one cell per image
% -------------------------------------------------------------------------
fprintf('\nMulti-threshold\n_____________________________________________\n');
fprintf('Retaining clusters contiguous with a significant cluster at p < %3.8f and k >= %3.0f\n\n', thresh(1), sizethresh(1));

for i = 1:length(thresh)
    
    dat = threshold(dat, thresh(i), 'unc', 'k', sizethresh(i), 'noverbose');
    dat = replace_empty(dat);
    
    for j = 1:nimgs
        
        sig{j}(:, i) = dat.sig(:, j);
        
    end
    
    
end

% prune - eliminate blobs not contiguous with significant blob at higher threshold

for j = 1:nimgs
    
    for i = 2:length(thresh)
        
        sig{j}(:, i) = iimg_cluster_prune(sig{j}(:, i), sig{j}(:, i-1), dat.volInfo);
        
    end
    
end


% display
for j = 1:nimgs
    
    fprintf('Image %3.0f of %3.0f\n_____________________________________________\n', j, nimgs);
    
    for i = 1:length(thresh)
        
        fprintf('Threshold p < %3.8f and k >= %3.0f\t%3.0f significant voxels', thresh(i), sizethresh(i), sum(sig{j}(:, i)));
        if i == 1
            fprintf(' defining blobs\n');
        else
            fprintf(' contiguous with a significant blob at a higher threshold\n');
        end
        
    end
    
    fprintf('\n');
    
end

fprintf('\n');






% create regions
% -------------------------------------------------------------------------

%r = cell(1, nimgs);
r = cell(1, length(thresh));

for i = 1:length(thresh)
    for j = 1 %:nimgs  % only handles first image in region.m anyway
        
        dat.sig(:, j) = sig{j}(:, i);
        r{i} = region(dat, 'noverbose');
        
    end
    
    [pcl{i}, ncl{i}] = posneg_separate(reparse_continguous(r{i}));
    
end

if nodisplay
    return
end


% set up fmridisplay
% -------------------------------------------------------------------------
useexisting = false;

if doPlot && ~isempty(o2) && isa(o2, 'fmridisplay')
        
        useexisting = true;
        
elseif doPlot && ~isempty(o2) && ~isa(o2, 'fmridisplay')
        error('o2 was passed in but is not an fmridisplay object.');

end

if doPlot && ~useexisting
    % Create one figure with multiple rows, for each image:
    
    o2 = canlab_results_fmridisplay([], 'multirow', nimgs);
    
end

% plot each image
% -------------------------------------------------------------------------

indx = 1;


for j = 1:nimgs
    
    if doPlot
        fprintf('\nMontage %d of %d\n_____________________________________\n', j, nimgs);
    end
    datplot = dat;
    datplot.dat = dat.dat(:, j);
    datplot.p = dat.p(:, j);
    if ~isempty(dat.ste), datplot.ste = dat.ste(:, j); end
    
    if doPlot
        if useexisting
            % Re-plot each result on all existing montage(s)
            % o2 = removeblobs(o2); do not remove because we may want to
            % add to selected montages or over other blobs
            
        else
            % Create a new montage set (sagg/ax) and plot this image on it only
            
%             figure;
%             o2 = canlab_results_fmridisplay([], 'multirow', nimgs);
%             
            %             o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 6, 'noverbose');
            %             axh = axes('Position', [0.05 0.4 .1 .5]);
            %             o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh, 'noverbose');
            
        end
    end
    
    for i = length(thresh):-1:1
        
        datplot.sig = sig{j}(:, i);
        
        if doWrite
            out=datplot;
            out.dat(out.sig==0)=0;
            out.fullpath=[pwd filesep 'Image' num2str(j) '_MultiThresh_' out.type '_P' num2str(thresh(i)) '_k' num2str(sizethresh(i)) '.nii'];
            write(out);
        end
        
        mycolors = [negcolors(i) negcolors(i) poscolors(i) poscolors(i)];
        
        if doPlot
            if useexisting
                
                o2 = addblobs(o2, region(datplot, 'noverbose'), 'splitcolor', mycolors, 'noverbose', varargin{:});
                
            else
                
                o2 = addblobs(o2, region(datplot, 'noverbose'), 'splitcolor', mycolors, 'noverbose', 'wh_montages', [indx indx+1]);
                
                %                 o2 = addblobs(o2, region(datplot, 'noverbose'), 'splitcolor', mycolors, 'wh_montages', indx);
                %                 o2 = addblobs(o2, region(datplot, 'noverbose'), 'splitcolor', mycolors, 'wh_montages', indx+1);
                
            end
        end
    end % thresholds
    
    indx = indx + 2;
    if doPlot
        drawnow, snapnow % for publish, etc.
    end
end

end % function
