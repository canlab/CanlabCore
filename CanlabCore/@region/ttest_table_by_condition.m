function [roi_table, region_obj] = ttest_table_by_condition(region_obj, DATA_OBJ, varargin)
% Extract region-by-region data from DATA_OBJ fmri_data object(s) in
% regions specified by region_obj, and return a table of means and stats,
% along with extracted data.
%
% First line: One-line summary description of function
%
% :Usage:
% ::
%
%     [roi_table, region_obj] = ttest_table_by_condition(region_obj, DATA_OBJ, [optional inputs])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2017 Tor Wager
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
%   **region_obj:**
%        A region object from CANlab core tools
%
%   **DATA_OBJ:**
%        an fmri_data object or cell vector of fmri_data objects with data
%        to extract.  Does not need to be in the same space as region
%        object -- extract_data.m method converts using mm coordinates.
%        Reslicing can induce variation in results due to interpolation.
%
% :Optional Inputs:
%   **'conditionnames', 'condition_names':**
%        Followed by cell vector of names for each condition
%        No spaces or special characters.
%
%
% :Outputs:
%
%   **roi_table:**
%        A table object containing ROI coordinates, names, and indices for
%        each region, along with means, and t-stats
%
%   **region_obj:**
%        The region object you entered in, with appended data,
%        images x conditions, stored in .dat field for each region.
%
% :Examples:
% ::
%
%    % give examples of code here
%    param1 = abc();
%    param2 = xyz();
%    [out1,out2] = func_call(param1, param2)
%
% :References:
%   CITATION(s) HERE
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

% Get volume from a region
volfunction = @(x) x.numVox .* prod(abs(diag(x.M(1:3, 1:3))));

if ~iscell(DATA_OBJ)
    DATA_OBJ = {DATA_OBJ};
end

nconditions = length(DATA_OBJ);
conditionnames = cell(1, nconditions);

for i = 1:length(DATA_OBJ)
    conditionnames{i} = sprintf('Cond_%d', i);
end

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case {'conditionnames', 'condition_names', 'conditions'}
                conditionnames = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if length(conditionnames) ~= nconditions
    error('Condition names and DATA_OBJ must be cells of same length');
end

nrois = length(region_obj);


% Initial table
% ------------------------------------------------------------------------
Index = (1:nrois)';
Name_of_Region = char(region_obj.shorttitle);
xyz_mm = cat(1, region_obj.mm_center);

for i = 1:nrois
    volume_mm3(i, 1) = volfunction(region_obj(i));
end

roi_table = table(Index, Name_of_Region, xyz_mm, volume_mm3);

% Initial table
% ------------------------------------------------------------------------

conditionnames = strrep(conditionnames, ' - ', '_vs_');
conditionnames = strrep(conditionnames, '-', '_vs_');
conditionnames = strrep(conditionnames, ' ', '_');

% extract ROI average data from each region
% Make regions x conditions matrix
% ------------------------------------------------------------------------

stars_by_condition = cell(nrois, nconditions);
[means_by_condition, t_by_condition, p_by_condition] = deal(zeros(nrois, nconditions));

for i = 1:nconditions
    
    rois_tmp = extract_data(region_obj, DATA_OBJ{i});
    
    % cat within regions to return data in region_obj object
    clear mydat
    for j = 1:nrois
        
        mydat = rois_tmp(j).dat(:, 1);
        
        % pad with longer of two
        [Anew Bnew] = padwithnan(region_obj(j).dat, mydat, 1);
        
        % region_obj(j).dat(:, i)
        region_obj(j).dat = [Anew Bnew];
        
    end
    
    % cat into image x roi matrix
    % for statistical tests
    dat = cat(2, rois_tmp.dat);
    
    % t-tests on each region across images (often subjects)
    [h, p, ci, stat] = ttest(dat);
    
    means_by_condition(:, i) = nanmean(dat)';
    
    t_by_condition(:, i) = stat.tstat';
    
    p_by_condition(:, i) = p';
    
    % Stars for each region
    for j = 1:length(p)
        
        if p(j) < .0015, mystr = '***';
        elseif p(j) < .015, mystr = '**';
        elseif p(j) < .055, mystr = '*';
        elseif p(j) < .105, mystr = '+';
        else mystr = ''; xadj = 0;
        end
        
        stars_by_condition{j, i} = mystr;
        
    end % loop through regions
    
    roi_table.(['t_' conditionnames{i}]) =  t_by_condition(:, i);
    
    signame = sprintf('sig%d', i);
    roi_table.(signame) = stars_by_condition(:, i);
    
end % loop through conditions

roi_table.means_by_condition = means_by_condition;
roi_table.p_by_condition = p_by_condition;


end % function
