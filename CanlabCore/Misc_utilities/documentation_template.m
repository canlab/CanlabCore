% Standard text for function documentation
%
% First line: One-line summary description of function
%
% Usage:
% -------------------------------------------------------------------------
% [list outputs here] = function_name(list inputs here, [optional inputs])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% Author and copyright information:
% -------------------------------------------------------------------------
%     Copyright (C) <year>  <name of author>
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
% xxx           xxx
%
% Outputs:
% -------------------------------------------------------------------------
% xxx           xxx
%
% Examples:
% -------------------------------------------------------------------------
%
% give examples here
%
% See also:
% * list other functions related to this one, and alternatives*

% Programmers' notes:
% List dates and changes here, and author of changes

% BELOW IS A STANDARD TEMPLATE FOR DEFINING VARIABLE (OPTIONAL) INPUT
% ARGUMENTS. MANY FUNCTIONS NEED TO PARSE OPTIONAL ARGS, SO THIS MAY BE
% USEFUL.


% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

% Defaults
% -----------------------------------
% initalize optional variables to default values here.
rowsz = [];
doplot = 0;
basistype = 'spm+disp';


% optional inputs with default values
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'rows', rowsz = varargin{i+1}; varargin{i+1} = [];
            case 'plot', doplot = 1; 
            case 'basistype', basistype = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% THIS IS ANOTHER WAY THAT IS MORE STREAMLINED, BUT ACTUALLY LESS INTUITIVE

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

% Defaults
% -----------------------------------
% Set color maps for + / - values
poscm = colormap_tor([1 0 .5], [1 1 0], [.9 .6 .1]);  %reddish-purple to orange to yellow
negcm = colormap_tor([0 0 1], [0 1 1], [.5 0 1]);  % cyan to purple to dark blue

% optional inputs with default values
% -----------------------------------
% - allowable_args is a cell array of argument names
% - avoid spaces, special characters, and names of existing functions
% - variables will be assigned based on these names
%   i.e., if you use an arg named 'cl', a variable called cl will be
%   created in the workspace

allowable_args = {'cl', 'ycut_mm', 'pos_colormap', 'neg_colormap', ...
    'surface_handles', 'existingfig'};

default_values = {[], [], poscm, negcm, ...
    [], 0};

% define actions for each input
% -----------------------------------
% - cell array with one cell for each allowable argument
% - these have special meanings in the code below
% - allowable actions for inputs in the code below are: 'assign_next_input' or 'flag_on'

actions = {'assign_next_input', 'assign_next_input', 'assign_next_input', 'assign_next_input', ...
    'assign_next_input', 'flag_on'};

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


% Below are some default headers

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Sub-functions
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



