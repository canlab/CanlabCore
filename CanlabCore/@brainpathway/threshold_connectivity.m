function threshold_connectivity(b_para, varargin)
% Threshold connnectivity matrices stored within a brainpathways object
%
% :Usage:
% ::
%
%     threshold_connectivity(b_para, [optional inputs])
%
% For objects: Type methods(object_name) for a list of object methods
% This is a handle object, so operations here modify the object in the
% workspace without passing it out.
%
%
% :Inputs:
%
%   **b_para:**
%        A brainpathway object
%
% :Optional Inputs:
%   **method:**
%        Correction method:
%        'fdr'  [default]
%        'unc'  uncorrected
%        'bonf' bonferroni based on num elements
%
%   **target:**
%       Which matrix to apply the threshold to:
%       'regions' [default]
%       'nodes'
%
%   **threshold:**
%       p or q-value threshold, numeric
%
% :Outputs:
%
%   b_para.connectivity fields are updated
%   e.g., b_para.connectivity.regions by default, unless you specify
%   'nodes'
% 
%   Arguments/threshold values are saved in:
%   b_para.connectivity.regions.r_thr_args
%   the equiv_p_threshold field is the P-value threshold after correcting
%   (FDR or bonf)
%
% :Examples:
% ::
%
%   Threshold at q < 0.05 FDR
%   threshold_connectivity(b_para)
%
%   Threshold at q < 0.01 FDR
%   threshold_connectivity(b_para, 'threshold', .01)
%
%   Threshold at p < 0.05 Bonferroni 
%   threshold_connectivity(b_para, 'method', 'bonf')
%
%   Threshold at p < 0.001 uncorrected 
%   threshold_connectivity(b_para, 'method', 'unc', 'threshold', .001)
%
%   Create a graph object from this and plot:
%   G = graph(b_para.connectivity.regions.r_thr, b_para.region_atlas.labels, 'upper', 'OmitSelfLoops');
%   figure; plot(G);
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2020 Tor Wager
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

% ----------------------------------------------------------------------
% Set defaults and parse inputs
% ----------------------------------------------------------------------
ARGS = parse_inputs(varargin{:});

% Logical flags
% ----------------------------------------------------------------------
if any(strcmp(varargin, 'noverbose')), doverbose = false; end
if any(strcmp(varargin, 'noplots')), doplots = false; end

if any(strcmp(varargin, 'fdr')), ARGS.method = 'fdr'; end
if any(strcmp(varargin, 'unc')), ARGS.method = 'unc'; end
if any(strcmp(varargin, 'bonf')), ARGS.method = 'bonf'; end

if any(strcmp(varargin, 'nodes')), ARGS.target = 'nodes'; end

% ----------------------------------------------------------------------
% Prep matrices
% ----------------------------------------------------------------------

r = b_para.connectivity.(ARGS.target).r;
p = b_para.connectivity.(ARGS.target).p;

trilp = tril(p, -1);
wh = logical(tril(ones(size(p)), -1)); % Select off-diagonal values (lower triangle)
trilp = double(trilp(wh));             % Vectorize and enforce double
trilp(trilp < 10*eps) = 10*eps;        % Avoid exactly zero values

% ----------------------------------------------------------------------
% Threshold depending on method
% ----------------------------------------------------------------------

switch ARGS.method
    
    case 'fdr'
% FDR-correct p-values across unique elements of the matrix

thr = FDR(trilp, ARGS.threshold);

if isempty(thr), thr = -Inf; end


    case 'unc'
        
        thr = ARGS.threshold;
        
    case 'bonf'
        
        k = size(r, 2);
        
        thr = ARGS.threshold ./ (k * (k - 1) ./ 2);
        
    otherwise error('threshold_connectivity: Unknnown threshold method')
        
end

        
sig = p < thr;


r_thr = r .* sig;

b_para.connectivity.(ARGS.target).r_thr = r_thr;

ARGS.equiv_p_threshold = thr;
b_para.connectivity.(ARGS.target).r_thr_args = ARGS;




end % Main function


function ARGS = parse_inputs(varargin)

p = inputParser;

% Validation functions - customized for each type of input
% ----------------------------------------------------------------------

valfcn_scalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar'});
valfcn_number = @(x) validateattributes(x, {'numeric'}, {'nonempty'}); % scalar or vector

% Validation: Region object, structure, or [x1 x2 x3] triplet
valfcn_custom = @(x) isstruct(x) || isa(x, 'region') || (~isempty(x) && all(size(x) - [1 3] == 0) && all(isnumeric(x)));

% Validation: [x1 x2 x3] triplet
valfcn_xyz = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'size', [1 3]});

valfcn_logical = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', '>=', 0, '<=', 1}); % could enter numeric 0,1 or logical
valfcn_text = @(x) validateattributes(x, {'char'}, {'nonempty'}); % text


% Required inputs
% ----------------------------------------------------------------------
% p.addRequired('x', valfcn_custom);
% p.addRequired('y', valfcn_custom);

% Optional inputs
% ----------------------------------------------------------------------
% Pattern: keyword, value, validation function handle

p.addParameter('method', 'fdr', valfcn_text);
p.addParameter('target', 'regions', valfcn_text);
p.addParameter('threshold', .05, valfcn_number);

% Parse inputs and distribute out to variable names in workspace
% ----------------------------------------------------------------------
% e.g., p.parse([30 1 0], [-40 0 10], 'bendpercent', .1);
p.parse(varargin{:});

ARGS = p.Results;

end % parse_inputs
