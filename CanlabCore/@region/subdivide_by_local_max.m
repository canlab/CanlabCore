function subregions = subdivide_by_local_max(r, varargin)
% Subdivide regions into sub-regions based on local peak Z-scores/maxima
%
% :Usage:
% ::
%
%    subregions = subdivide_by_local_max(r, ['mm_distance', value], ['noorthviews'])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
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
% :Optional Inputs:
%
%   **mm_distance:**
%        Followed by mm distance minimum for dividing subclusters
%
%   **noorthviews:**
%        Suppress display of orthviews
%
% :Outputs:
%
%   **subregions:**
%        subdivided region object
%
% :See also:
% region.subdivide_by_atlas, subclusters_from_local_max, cluster_local_maxima
%
% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%
%    BELOW IS A STANDARD TEMPLATE FOR DEFINING VARIABLE (OPTIONAL) INPUT
%    ARGUMENTS. MANY FUNCTIONS NEED TO PARSE OPTIONAL ARGS, SO THIS MAY BE
%    USEFUL.
% ..

% ..
%    DEFAULTS AND INPUTS
% ..
subregions = region();

% optional inputs with default values
% -----------------------------------
% - allowable_args is a cell array of argument names
% - avoid spaces, special characters, and names of existing functions
% - variables will be assigned based on these names
%   i.e., if you use an arg named 'cl', a variable called cl will be
%   created in the workspace

allowable_args = {'mm_distance' 'noorthviews'};

default_values = {10, 0};

% define actions for each input
% -----------------------------------
% - cell array with one cell for each allowable argument
% - these have special meanings in the code below
% - allowable actions for inputs in the code below are: 'assign_next_input' or 'flag_on'

actions = {'assign_next_input', 'flag_on'};

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


for i = 1:length(r)
    subcl{i} = subclusters_from_local_max(r(i), mm_distance);
end

subcl = cat(2, subcl{:});

for i = 1:length(subcl)
    subregions(i) = cluster2region(subcl(i));
    
    subregions(i).val = subregions(i).Z';
    subregions(i).dim = r(1).dim;
    
end

if ~noorthviews
    orthviews(subregions, 'unique');
end

end % function

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Sub-functions
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



