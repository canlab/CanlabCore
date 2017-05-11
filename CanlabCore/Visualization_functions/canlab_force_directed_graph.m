function [stats, handles] = canlab_force_directed_graph(activationdata, varargin)
% Creates a force-directed graph from a set of variables, and plots
% clusters on 3-D brain as well if entered. Requires matlab BGL toolbox.
%
% :Usage:
% ::
%
%    canlab_force_directed_graph(activationdata OR connection matrix, ['cl', cl])
%
%
% ..
%     Author and copyright information:
%     -------------------------------------------------------------------------
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
%   **activationdata:**
%        observations x variables matrix of data to be
%        inter-correlated
%
%        OR
%
%        signed, thresholded connection matrix to be used (e.g.,
%        thresholded t-values from multi-subject group analysis
%
% :Optional Inputs: Enter keyword followed by variable with values
%
%   **'cl':**
%        followed by clusters or region structure with brain clusters
%
%   **'threshtype':**
%        followed by threshold type; 'bonf' is option now
%
%   **'connectmetric':**
%        followed by node connection metric 'corr' or 'partial_corr'
%
%   **'sizescale':**
%        Followed by values to use in sizing of nodes on graph
%        'linear' 'sigmoidal' [default] or 'custom'
%
%   ** 'sizes' **
%       Followed by vector of point sizes for all points - works ONLY if sizescale
%       is set to 'custom'
%
%   **'setcolors':**
%        Cell array of colors for each group, [1 x g]
%
%   **'rset':**
%        Cell of vectors, with indices (integers) of member
%        elements in each group, [1 x g] cell
%        rset can ALSO be a vector of integers, i.e., output
%        from clusterdata
%
%   **'names':**
%       followed by cell array of names for each region/object
%
%   **'namesfield':**
%       followed by name of field to extract region/object names from in
%       clusters structure
%
%   **'linewidth':**
%       followed by line width
%
%   **'linestyle':**
%       followed by line style. Default is 'curved', anything else is
%       'straight'
%
% :Output:
%
%   **stats:**
%        structure with descriptive statistics, including
%        betweenness-centrality, degree of each node
%
% :Examples:
% ::
%
%    [stats, handles] = canlab_force_directed_graph(activationdata, 'cl', cl, 'namesfield', 'shorttitle');
%    [stats, handles] = canlab_force_directed_graph(activationdata, 'cl', cl, 'namesfield', 'shorttitle', 'degree');
%    [stats, handles] = canlab_force_directed_graph(activationdata, 'cl', cl, 'namesfield', 'shorttitle', 'degree', 'rset', rset, 'setcolors', setcolors);
%
%    [stats, handles] = canlab_force_directed_graph(G, 'sizescale', 'custom', 'sizes', [6*ones(17, 1); 12*ones(7, 1)]);

% ..
%    DEFAULTS AND INPUTS
% ..

shan = [];              % outputs
spherehan = [];
handles = [];

cl = [];
threshtype = 'bonf';
connectmetric = 'corr';  % or partial_corr
sizescale = 'sigmoid';
sizes = ones(size(activationdata, 1));
linestyle = 'curved';

setcolors = [];         % control of color subgroups
rset = [];

ptsizetype = 'bc';
names = [];
namesfield = [];        % enter field name, e.g., 'shorttitle'
linewidth = 1;

% Variable arg inputs
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            % reserved keywords
            case 'degree', ptsizetype = 'degree';
                %case 'design'
                
                % functional commands
            case 'cl', cl = varargin{i + 1}; varargin{i + 1} = [];
                
            case 'linewidth', linewidth = varargin{i + 1}; varargin{i + 1} = [];
                
            case 'linestyle', linestyle = varargin{i + 1}; varargin{i + 1} = [];
                
            case {'threshtype' 'connectmetric' 'sizescale' 'setcolors' 'rset' 'names' 'namesfield' 'sizes'}
                eval([varargin{i} ' = varargin{i + 1}; varargin{i + 1} = [];'])
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Check vars, etc.
% ---------------------------------------------------------------------
if isdir('/Users/tor/Documents/matlab_code_external/matlab_bgl')
    g = genpath('/Users/tor/Documents/matlab_code_external/matlab_bgl');
    addpath(g);
end

if ~exist('fruchterman_reingold_force_directed_layout.m', 'file')
    error('Must have Matlab BGL toolbox on path (external toolbox)');
end

activationdata = double(activationdata);

if ~isempty(cl) && (size(activationdata, 2) ~= length(cl))
    disp('activationdata must have as many columns as clusters has elements.')
    error('check activationdata and clusters to make sure they match.')
end

if isempty(rset)
    rset = {1:size(activationdata)};
end

if isempty(setcolors)
    setcolors = scn_standard_colors(length(rset));
end

if ~isempty(namesfield)
    names = {cl.(namesfield)}';
end

if ~iscell(rset) % then integer vector
    for i = 1:max(rset)
        rset2{i} = find(rset == i);
    end
    rset = rset2;
end


% -------------------------------------------------------------------------
% CALCULATIONS
% -------------------------------------------------------------------------
if issymmetric(activationdata)
    fprintf('Thresholded association matrix detected.\n')
    [b, C, r, sig] = deal(activationdata);
    thr = NaN;
    
else
    fprintf('activationdata appears to be raw data. Running inter-correlations. See also xcorr_multisubject.\n')
    
    [r, rp] = corr(double(activationdata));
    sz = size(r, 1);
    
    switch threshtype
        case 'bonf'
            thr = .05 ./ (sz * (sz - 1) / 2);  % bonferroni...
    end
    
    % Partial correlations
    [b, p] = calc_partial_r(activationdata);
    
    % Threshold
    sig = sparse(double(p < thr & b > 0));
    signeg = sparse(-double(p < thr & b < 0));
    sig = sig + signeg;
    
    
    % Threshold based on significant partial regression effects
    % C is connectivity matrix for graph
    
    switch connectmetric
        case 'partial_corr'
            
            C = b;
            C(~sig) = 0;
            
        case 'corr' % default
            C = r;
            C(~sig) = 0;
            
        otherwise error('Unknown connectmetric');
    end
    
    % end if issymmetric
end

% Enforce format for graph
C = (C + C') ./ 2;
C = sparse(C);

% Graph Stats
% -----------------------------------------------
bc = betweenness_centrality(abs(C));  % abs if ignoring neg connections

% shortest paths: used as estimate of connectivity
[D S] = mean_path_length(r, rset);

deg = full(sum(C ~= 0))';

stats.names = names;
stats.rset = rset;
stats.C = C;
stats.r = r;
stats.b = b;
stats.thr = thr;
stats.betweenness = bc;
stats.path_length_distance = D;
stats.degree = deg;
stats.mean_path_by_rset = S;


switch ptsizetype
    case 'bc'
        ptsizemetric = bc;
    case 'degree'
        ptsizemetric = deg;
    case 'none'
        ptsizemetric = ones(size(bc));
    otherwise error('Unknown ptsizetype');
end
if strcmp(sizescale,'custom')
   fprintf('Node size reflects custom scaling\n');
else
    fprintf('Node size reflects %s\n', ptsizetype);
end
% -------------------------------------------------------------------------
% Force-directed graph
% -------------------------------------------------------------------------

% Xc is coordinate matrix
Xc = fruchterman_reingold_force_directed_layout(C);

if ~isempty(cl)
    create_figure('graph', 1, 2);
else
    create_figure('graph');
end

switch lower(linestyle)
    
    case 'curved'
        
        han = nmdsfig_tools('drawlines',Xc, sig, [.5 .5 .5; 0 1 1],[{'-'} {':'}], .1);
        handles.lh_pos = han.hhp;
        handles.lh_neg = han.hhn;% Dots
        %set(handles.lh_pos, 'Color', [.5 .5 .5])
        set(handles.lh_neg, 'Color', [0 .5 1], 'LineStyle', ':')

    otherwise
        
        gplot(sig>0, Xc, 'k');
        axis image
        lh = findobj(gca, 'Type', 'line');
        set(lh, 'Color', [.5 .5 .5])
        
        handles.lh_pos = lh;
        
        gplot(sig<0, Xc, 'b');
        lh = findobj(gca, 'Type', 'line', 'color', 'b');
        set(lh, 'Color', [0 .5 1], 'LineStyle', ':')
        
        handles.lh_neg = lh;
        
end

switch sizescale
    case 'linear'
        %         sizes = 4 + 15 * intensity.(cnames{i}) ./max(intensity.(cnames{i}));  % LC none of these variables exist in the workspace
        sizes = 4 + 15 * ptsizemetric ./ max(ptsizemetric);
        
    case 'sigmoid'
        % Rationale: avoids some VERY large points in areas highly
        % focused on as a priori ROIs, with extreme z-scores in
        % intensity relative to other areas.
        
        sizes = sigmoidscale(ptsizemetric);
        
    case 'custom'
        % we already have sizes
end

ph = [];
texth = {};

for k = 1:length(rset)
    
    texth{k} = [];
    
    for j = 1:length(rset{k})
        
        ph(end+1) = plot(Xc(rset{k}(j), 1), Xc(rset{k}(j), 2), 'o', 'Color', setcolors{k} ./ 2, 'MarkerSize', sizes(rset{k}(j)), 'MarkerFaceColor', setcolors{k}, 'LineWidth', linewidth);
        
        if ~isempty(names)
            offset = .02 * range(get(gca, 'XLim'));
            texth{k}(end+1) = text(Xc(rset{k}(j), 1)+offset, Xc(rset{k}(j), 2), strrep(names{rset{k}(j)}, '_', ' '), 'Color', 'k', 'FontSize', 14);
        end
    end
    
end

handles.ph = ph;
handles.texth = texth;

xlim = get(gca, 'XLim');  rg = range(xlim) * [-.05 .05]; xlim = xlim + rg;
ylim = get(gca, 'YLim'); rg = range(ylim) * [-.05 .05]; ylim = ylim + rg;
set(gca, 'XLim', xlim, 'YLim', ylim);
axis off
drawnow

stats.Xc = Xc;
stats.sizes = sizes;




if isempty(cl), return, end

% -------------------------------------------------------------------------
% Subcortical surface
% Only if cl is entered
% -------------------------------------------------------------------------

% NOTE: you need 3dHeadUtilityLite on your matlab path.
subplot(1, 2, 2)

% regioncenters
xyz = cat(1, cl.mm_center);

DB = struct('xyz', xyz, 'x', xyz(:, 1), 'y', xyz(:, 2), 'z', xyz(:, 3));

% sizes = sigmoidscale(ptsizemetric, 2, 6);
sizes = sizes ./ 2.5; % rescale to avoid huge ones

% Make Brain with Spheres
% ------------------------------------------------------------

[shan, spherehan] = connectivity3dbrain(xyz, rset, sizes, setcolors, names);

% Add lines
% ------------------------------------------------------------

[~, linehandles] = cluster_nmdsfig_glassbrain(cl,ones(length(cl), 1), {[.5 .5 .5]}, sig, [], 'samefig', 'nobrain', 'noblobs', linestyle);
set(linehandles, 'LineWidth', 2)

drawnow
%saveas(gcf, fullfile(savefigdir, ['graph3d_v2_' cnames{i} '.png']));

handles.surfhan = shan;
handles.spherehan = spherehan;
handles.linehandles = linehandles;

end % function




function [b p] = calc_partial_r(activationdata)

X = activationdata;
X = zscore(X);
%X(isnan(X)) = mean(X(~is;

clear b p

for i = 1:size(X, 2)
    y = X(:, i);
    xx = X;
    xx(:, i) = 1;  % intercept; need it, and also placeholder
    
    [bb, dev, statsglm] = glmfit(xx, y, 'normal', 'constant', 'off');
    
    b(:, i) = bb;
    p(:, i) = statsglm.p;
    p(i, i) = 1;
    b(i, i) = NaN;
    
end


end % function





function sizes = sigmoidscale(sizes, varargin)
% sizes = sigmoidscale(sizes, [lower bound], [upper bound])
%
% rescales a vector based on sigmoid function of zcore(input values)

% Rationale: avoids some VERY large points in areas highly
% focused on as a priori ROIs, with extreme z-scores in
% intensity relative to other areas.

% scale size of nodes - intensity, or betweenness
sizes = zscore(sizes);

A = 4; % lower asymptote
K = 15; % upper asymptote

if nargin > 2
    A = varargin{1};
    K = varargin{2};
end

B = 2;   % growth rate
v = .5;  % high-growth asymptote
Q = .5;
M = .5;  % time of max growth, if v = Q

richards = @(x, A, K, B, v, Q, M) A + (K - A) ./ ((1+Q*exp(-B*(x-M))).^(1/v));

%figure; plot(sort(sizes), richards(sort(sizes), A,K,B,v,Q,M));
sizes = richards(sizes, A,K,B,v,Q,M);

end




function [C S] = mean_path_length(r, rset)
% compute matrix of 1/path lengths for n x n matrix of regions, C (connectivity)
% and 1 / mean path length for sets of regions specified by rset
%
% r is a correlation matrix
% rset is a cell array of length k, for k sets, with vectors describing the
% indices of members of each set.
%
% e.g., r = region_r{i};

r(isnan(r)) = 0;
r = (r' + r) ./ 2;  % enforce symmetry, just in case


% shortest paths: used to calculate connectivity
rtmp = sparse(r);
rtmp(rtmp < 0) = 0;

% Note: Uses Floyd-Warshall method if 10% non-zero elements or more
C = all_shortest_paths(rtmp);


% for each set, average elements of D corresponding to pairs of SETS
% return averages in S
%
% e.g., if D is min path length, S is the average min path length for Set 1
% vs. 2, 1 vs. 3, 2 vs. 3, etc.

% for inf values, we must impute some finite value or every average is Inf
% (i.e., if any regions are unconnected)
% impute max:
% but maybe better to work on unthresholded r matrix.
% Or, even better, return 1/S, 1 / mean path length, so unconnected
% regions/sets get 0.
%  mx = max(D(~isinf(D)));
%  D(isinf(D)) = mx;

S = zeros(length(rset));

for m = 1:length(rset)-1
    for n = m:length(rset)
        
        % set m to set n relationships (e.g., shortest path lengths)
        vals = 1 ./ C(rset{m}, rset{n});
        
        % average, excluding "self-connections" with value Inf
        S(m, n) = mean(vals(~isinf(vals)));
        
    end
end

S = S + S';

end



function [shan, spherehan] = connectivity3dbrain(regioncenters, rset, sizes, setcolors, names)

% Colors
% ----------------------------------------------------------

ctxcolors = cell(1, size(regioncenters, 1));

for i = 1:length(rset)
    
    for j = 1:length(rset{i})
        
        wh = rset{i}(j);
        
        ctxcolors(rset{i}) = setcolors(i);
        
    end
end

% Brain
% ----------------------------------------------------------

%create_figure('subcortex');

shan = addbrain('limbic');
delete(shan(end));
shan = shan(1:end-1);

shan = [shan addbrain('hires left')];
shan = [shan addbrain('brainstem')];

set(shan(end-1), 'FaceColor', [.5 .5 .5], 'FaceAlpha', .15);
set(shan(1:end-2), 'FaceColor', [.5 .5 .5], 'FaceAlpha', .3);
set(shan(end), 'FaceColor', [.5 .5 .5], 'FaceAlpha', .3);

view(89, 1);
lightRestoreSingle;

lighting gouraud;

% Spheres
% ----------------------------------------------------------

% add spheres
% cortex is special, because it involves different colors

% xyz = regioncenters(rset{1}, :);
% sz = sizes(rset{1});
% if ~isempty(names)
%     mynames = names(rset{1});
% end

xyz = regioncenters;
sz = sizes;
if ~isempty(names)
    mynames = names;
end

for i = 1:size(xyz, 1)
    spherehan(i) = cluster_image_sphere(xyz(i, :), 'color', ctxcolors{i}, 'radius', sz(i));
    
    if ~isempty(names)
        offset = .04 * range(get(gca, 'XLim'));
        text(xyz(i, 1)+offset, xyz(i, 2)+offset, xyz(i, 3)+offset, mynames{i}, 'Color', 'k', 'FontSize', 14);
    end
end
spherehan = {spherehan};

% for i = 2:length(rset)
%
%     xyz = regioncenters(rset{i}, :);
%     sz = sizes(rset{i});
%     if ~isempty(names)
%     mynames = names(rset{i});
%     end
%
%     spherehan{i} = cluster_image_sphere(xyz, 'color', setcolors{i}, 'radius', sz);
%
%     if ~isempty(names)
%         offset = .02 * range(get(gca, 'XLim'));
%         for j = 1:size(xyz, 1)
%             text(xyz(j, 1)+offset, xyz(j, 2)+offset, xyz(j, 3)+offset, mynames{j}, 'Color', 'k', 'FontSize', 14);
%         end
%     end
%
% end


end % function
