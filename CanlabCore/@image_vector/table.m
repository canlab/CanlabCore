function [region_obj, results_table, varargout] = table(w, varargin)
% Print a table of all regions in a thresholded image. Return a Matlab table object
% and labeled region object. Used with fmri_data, statistic_image, and
% image_vector objects.
%
% This function has two "modes":
% - Without an optional atlas object input, it uses region.table to generate a table based
% on contiguous clusters (separated by those with positive and negative peak values) 
%
% - With an optional atlas object input, it uses image_vector.subdivide_by_atlas to 
% generate a table based on atlas regions covered by the input image. The table has
% one row per atlas region. The region object is also divided by the atlas,
% so one large contiguous blob in the image may cover multiple rows with
% different labeled regions.
%
% :Usage:
% ::
%
%    [poscl, negcl] = table(cl, [optional inputs])
%
% - By default, region.table() separates clusters into subregions with positive
%   and negative values. Thus, the number of rows may not match the original number of regions. 
%   To turn this feature off, use 'nosep'.
% - By default, region.table() re-sorts the regions to group the table rows by macro-scale brain
%   structures (cortex, basal ganglia, etc.). So the ordering in the table may not match
%   the original region object. To turn this feature off, use 'nosort'.
%
% :Optional inputs:
%
%   **atlas_obj**
%       Any atlas-class object, e.g., loaded by load_atlas().
%
%   **k:**
%        Print only regions with k or more contiguous voxels
%
%   **nosep:**
%        do not separate cl with pos and neg effects based on peak in .val
%
%   **names:**
%        name clusters manually before printing to table and output; saves in .shorttitle field
%
%   **forcenames:**
%        force manual naming of cl by removing existing names in .shorttitle field
%
%   **nosort:**
%        Do not sort rows by network/brain lobe [default is to sort]
%
%   **legacy:**
%        force manual naming of cl by removing existing names in .shorttitle field
%
%   **nolegend:**
%        omit table legend
%
% :Outputs:
%
%   Returns region objects for cl with pos and neg effects
%   - autolabeled if Neuroimaging_Pattern_Masks and atlas tools are available on Matlab path
%   - limited by size if entered
%   - manually named if requested (see optional inputs)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2011  Tor Wager
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
% :Examples:
% -------------------------------------------------------------------------
% ::
% Example 1:
% % Complete group analysis of a standard dataset
% % Do analysis and prep results region object:
%
%   img_obj = load_image_set('emotionreg');         % Load a dataset
%   t = ttest(img_obj, .005, 'unc');                % Do a group t-test
%   t = threshold(t, .005, 'unc', 'k', 10);         % Re-threshold with extent threshold of 10 contiguous voxels
%   r = region(t);                                  % Turn t-map into a region object with one element per contig region
%
%   Label regions and print a table:
%   [r, region_table, table_legend_text] = autolabel_regions_using_atlas(r);  
%                                                   % Label regions. Can be skipped because 'table' below attempts to do this automatically
%   table(r);                                       % Print a table of results using new region names
%
%   [rpos, rneg] = table(r);                        % Print and table and return region object separated into regions with positive vs. negative statistic values (from .Z field)

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..
%    July 2018:  Autolabel update and "new 2018 version", Tor Wager. Also added legend text.

n_cols = 140;                       % 140 good for HTML reports
sep_str = repmat('_', 1, n_cols);   % see textwrap

k = 0;
dosep = true;
donames = false;        % name clusters before printing to table and output; saves in .shorttitle field (legacy only)
forcenames = false;     % force naming of cl by removing existing names in .shorttitle field (legacy only)
dolegacy = false;
dosortrows = true;          % sort rows by area
dolegend = true;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            % functional commands
            case {'k', 'maxsize'}, k = varargin{i+1};
            case 'nosep', dosep = 0;
            case {'names', 'name', 'donames'}, donames = true;
            case 'forcenames', forcenames = true;
                
            case {'nosort', 'nosortrows'}, dosortrows = false;
                
            case {'legacy', 'dolegacy'}, dolegacy = true;
                
            case 'nolegend', dolegend = false;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

 

cl = region(w); % Convert w to cl; MS bugfix 09172024

if exist()

    atl = load_atlas('canlab2024');
end

[poscl, negcl, results_table] = table(cl, varargin)


[w_atlas, w_region_obj] = subdivide_by_atlas(w, atl);

% Make a table of labeled regions with significant voxels in the weight map

[region_obj, region_table] = autolabel_regions_using_atlas(w_region_obj, atl);
region_table = region_table(:, 1:3);
