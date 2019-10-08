function [stats, hh, hhfill, table_group, multcomp_group] = image_similarity_plot(obj, varargin)
% Associations between images in object and set of 'spatial basis function' images (e.g., 'signatures' or pre-defined maps)
% - Similarity metric is point-biserial correlations or cosine-similarity
%
% Usage:
% ::
%
%    stats = image_similarity_plot(obj, 'average');
%
% This is a method for an image_vector object that compares the spatial
% similarity of input image(s) to a specified set of a priori spatial basis maps.
% It returns similarity values (Pearson's r) to each a priori basis map,
% and a plot that shows these values.  If multiple images are entered, the
% 'average' function can return a plot with standard error bars and
% statistics on the significance of the correlation with each basis map
% (across input images) and the differences (inhomogeneity) in similarity across basis
% maps.  If a grouping variable is entered, statistics are calculated for
% multivariate differences across the groups of input images. Such
% differences are assessed treating the basis maps as variables and input
% images as cases.  The basis sets are "NPSplus" (the default), which
% includes the NPS map from Wager et al. 2013, Romantic Rejection
% classifier (Woo 2015), Negative emotion map (Chang 2015), and vicarious
% pain (Krishnan et al. 2016).  Other sets are "bucknerlab" including 7 cortical [only]
% networks from the Buckner Lab's 1000-person resting-state analyses and
% "kragelemotion", including 7 emotion-predictive maps from Kragel 2015.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2015 Tor Wager, stats added by Phil Kragel
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
%   **mapset**
%       Followed by one of the keywords below, or by an fmri_data object
%       containing maps you want to apply to (compare similarity with) input image objects
%       If you enter a custom image object, also enter 'networknames'
%       followed by cell array of names for each images in mapset, OR enter
%       names in mask image obj.image_names (these are used by default)
%
%   **bucknerlab**
%        Use 7 network parcellation from Yeo et al. as basis for comparisons
%        Cortex only.  BUT also:
%        'bucknerlab_wholebrain': 7 networks in cortex, BG, cerebellum
%        'bucknerlab_wholebrain_plus': 7 networks in cortex, BG, cerebellum
%        + SPM Anatomy Toolbox regions + brainstem
%
%   **kragelemotion**
%        Use 7 emotion-predictive models from Kragel & LaBar 2015 for
%        basis of comparisons
%
%   **allengenetics**
%        Five maps from the Allen Brain Project human gene expression maps
%        from Luke Chang (unpublished)
%
%   **bgloops**
%        5-basal ganglia parcels and 5 associated cortical
%        networks from Pauli et al. 2016.  Also 'pauli'
%        'bgloops17', 'pauli17' : 17-parcel striatal regions only from Pauli et al. 2016
%
%   **fibromyalgia**
%        3 neural classifiers used to predict FM in Lopez-Sola et al 2017
%        also 'fm','fibro'
%
%   **pain_pdm**
%          11 high-dimensional pain mediator maps. PDM1-10 are individual, 
%          orthogonal maps. The Joint PDM is a weighted combination of the 
%          other 10 PDMs.
% 
%
%   **End of mapset options**
%   
% 	**compareGroups**
%        Perform multiple one-way ANOVAs with group as a factor (one for
%        each spatial basis); requires group as subsequent input
%
%   **group**
%        Indicates group membership for each image
%
%   **plotstyle:**
%            'wedge' [default] or 'polar'
%
%   **noplot**
%        Omits plot (print stats only)
%
%   **nofigure**
%       Omit creation of new figure
%
%   **notable**
%       Omit printing of similarity table
%
%   **cosine_similarity**
%        Use cosine similarity instead of Pearson's r
%        Pearson's r is the default similarity metric.
%
%   **colors**
%             followed by cell array of colors, one for each image/image
%             group
%               Default colors for multi-line plots
%                       Color
%               1		Red
%               2		Green
%               3		Dark Blue
%               4		Yellow
%               5		Pink
%               6		Turquoise
%
%   **bicolor**
%        For wedge plot only, plot positive entries in groupColors{1} and
%        negative entries in groupColors{2}. Otherwise, will plot negative
%        entries with stripes on wedges [default].
%
%   **'dofixrange':**
%        Set min and max of circles numbers (values on polar axis)
%        Follow by range vector: [min_val max_val] OR
%        one radius value for wedge
%
%
% :Outputs:
%
%   **stats:**
%        Structure including:
%           - .r, Similarity (correlations or cosine sim) in [mask images x images in obj] matrix
%           - .t, T-test (if 'average' is specified)
%           - .line_handles Handles to polar plot lines so you can
%             customize
%           - .fill_handles Handles to polar plot fills so you can
%             customize
%           - .table_spatial, ANOVA table with subject as row factor and
%             spatial basis as column factor (one way repeated measures
%             ANOVA, requires 'average' to be specified)
%           - .multcomp_spatial, multiple comparisons of means across
%             different spatial bases, critical value determined
%             by Tukey-Kramer method (see multcompare)
%   **hh:**
%             Handles to lines
%
%   **hhfill:**
%             Handles to fill areas
%
%   **table_group**
%             multiple one-way ANOVA tables (one for each
%             spatial basis) with group as column factor (requires
%             'average' to be specified)
%
%   **multcomp_group**
%             mutiple comparisons of means across groups, one output
%             cell for each spatial basis, critical value determined
%             by Tukey-Kramer method (see multcompare)
%
%
% :Examples:
% ::
%
% testimgs = load_image_set('emotionreg');
% stats = image_similarity_plot(testimgs, 'cosine_similarity', 'bucknerlab');
% stats = image_similarity_plot(testimgs, 'cosine_similarity', 'bucknerlab', 'polar');
% stats = image_similarity_plot(testimgs, 'cosine_similarity', 'bucknerlab', 'plotstyle', 'polar', 'average');
% stats = image_similarity_plot(testimgs, 'cosine_similarity', 'bucknerlab', 'plotstyle', 'polar', 'average', 'colors', {[1 .5 0]});
% 
%    % corrdat is an fmri_data object with 18 images from searchlight
%    % correlation in it.  Then:
%    stats = image_similarity_plot_bucknermaps(corrdat, 'average');
%
%    % t_diff is a thresholded statistic_image object
%    stats = image_similarity_plot_bucknermaps(t_diff);
%
%    stats = image_similarity_plot(fmri_data(img), 'cosine_similarity', 'bucknerlab', 'colors', color);
%
% :See also:
%
% tor_polar_plot, tor_wedge_plot

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% List dates and changes here, and author of changes
% 11/30/2015 (Phil Kragel)
%   -   added anova (rm) comparing means across spatial bases
%   -   added anova (1-way) comparing means across groups for each spatial
%       basis (e.g., for each buckner network)
% 12/15/2015 (Phil Kragel)
%   - added option to omit plotting
% 5/10/2016 (Phil Kragel)
%   - added option to use cosine similarity instead of Pearson
% 8/21/2017 (Stephan Geuter)
%   - fixed header for printing similarity table
% 2017/09/07 Stephan Geuter
%   - added option for percent overlap of binary masks (see also
%   canlab_pattern_similarity.m and riverplot.m
%   Changed metric selection to string format.
% 2018/1/9  tor: changed default colors for compat with wedge plot,
% debugged wedge plot with average option.
% 2018/1/16 Stephan: added pain PDM mediators as mapsets
% 
% ..


% PRELIMINARIES
% ------------------------------------------------------------------------


% If mask is an atlas object, convert to fmri_data object containing the
% probability_map data

if isa(obj, 'atlas')
    obj = atlas_get_probability_maps(obj);
end

n_obs = size(obj.dat, 2); % number of images to test + plot


% DEFAULTS AND INPUTS
% ------------------------------------------------------------------------

[hh, hhfill] = deal(' ');
doaverage = false;       % initalize optional variables to default values here.
force_noaverage = false; % averaging mode determined by plot style below, which is problematic for some
                         % functions, e.g., riverplot

mapset = 'bucknerlab';  % 'bucknerlab'
table_group = {}; %initialize output
multcomp_group = {}; %initialize output
dofigure = true;
noplot = false;

groupColors = [{[1 .9 0] [0 0 1]}];  % for bicolor wedge defaults: yellow pos, blue neg  scn_standard_colors(n_obs);
morecolors = scn_standard_colors(n_obs + 2);
groupColors = [groupColors morecolors(5:end)];  % avoid redundancy

dofixRange = 0;
printTable = true;
% changed metric selection to string format (SG 2017/09/07)
sim_metric = 'corr'; % default: correlation
% doCorr = 1;
% doCosine = 0; %do cosine similarity
plotstyle = 'wedge'; % or 'polar'
bicolor = false;

% optional inputs with default values
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'average'
                doaverage = true;
                bicolor = true;  % if wedge, bicolor only, no lines.
               
            case 'noaverage'
                force_noaverage = true;
                
            case 'cosine_similarity', sim_metric = 'cosine';
                
            case 'binary_overlap', sim_metric = 'overlap';
                
            case {'bucknerlab', 'bucknerlab_wholebrain' 'bucknerlab_wholebrain_plus' ...
                    'kragelemotion' 'allengenetics' ...
                    'pauli' 'bgloops' 'pauli17' 'bgloops17' 'fm' 'fibro' 'fibromyalgia' ...
                    'bgloops_cortex' 'painsig' 'pain_pdm' 'pdm' }
                
                mapset = varargin{i};
                
            case 'mapset'
                %mapset = 'custom';
                mapset = varargin{i + 1}; varargin{i + 1} = [];
                
                %case 'basistype', basistype = varargin{i+1}; varargin{i+1} = [];
            case 'compareGroups'
                compareGroups = true;
                group = varargin{i+1};
                
            case 'noplot'; noplot = true;
                
            case 'nofigure', dofigure = false;
                
            case {'fixedrange', 'dofixrange'}
                dofixRange = 1;
                fixedrange = varargin{i+1};
                
            case 'colors'
                groupColors = varargin{i + 1}; varargin{i + 1} = [];
                
            case 'bicolor'
                bicolor = true;
                
            case 'networknames' % do nothing, handle later
                
            case 'notable'
                printTable = false;
                
            case 'plotstyle'
                plotstyle = varargin{i + 1}; varargin{i + 1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% ADJUSTMENTS to default behavior
% ------------------------------------------------------------------------

% Wedge handles inputs a bit differently from 'polar'.  The
% default input is regions x n observations to make a
% series of n polar line plots. but here we want averages,
% because it doesn't make sense to have a series of wedge
% plots. So we transpose r and average across observations
% by default.

% This change in default behavior is problematic for some applications,
% e.g., riverplots. 
if n_obs > 1 && strcmp(plotstyle, 'wedge') && ~force_noaverage
    
    % We have a wedge plot with multiple obs/images. Default to 'average'
    % mode, not multi-line-plot. For polar plots, do multi line plot.
    
    doaverage = 1;
    bicolor = true;  % if wedge, bicolor only, no lines.
end

% Colors look pretty funny with 'average' and 'polar', so adjust defaults
% if colors are not entered

if doaverage && strcmp(plotstyle, 'polar') && ~any(strcmp(varargin, 'colors'))
   
    groupColors = {[1 0 0]};  % unicolor. can change line and err with 2nd color entry.
    
end

% Load image set: Most of the hard work
% ------------------------------------------------------------------------

[mask, networknames, imagenames] = load_image_set(mapset); % Patterns/regions to apply.

% Re-load names if entered.
% ------------------------------------------------------------------------
wh = strcmp(varargin, 'networknames');
if any(wh)
    wh_names = find(wh) + 1;
    networknames = varargin{wh_names(1)};
    
    % Must be row vector of cells
    if iscolumn(networknames), networknames = networknames'; end
end



% Deal with space and empty voxels so they line up
% ------------------------------------------------------------------------

mask = replace_empty(mask); % add zeros back in

mask = resample_space(mask, obj);

obj = remove_empty(obj);
nonemptydat = ~obj.removed_voxels; % remove these

obj = replace_empty(obj);

% Correlation or other similarity metric
% ------------------------------------------------------------------------

% Point-biserial correlation is same as Pearson's r.
% Gene V. Glass and Kenneth D. Hopkins (1995). Statistical Methods in Education and Psychology (3rd edition ed.). Allyn & Bacon. ISBN 0-205-14212-5.
% Linacre, John (2008). "The Expected Value of a Point-Biserial (or Similar) Correlation". Rasch Measurement Transactions 22 (1): 1154.
% http://www.andrews.edu/~calkins/math/edrm611/edrm13.htm#POINTB

% If both binomial, could use Dice coeff:
%dice_coeff = dice_coeff_image(mask);

% if map or series of maps, point-biserial is better.

% If mask is an atlas object, convert to fmri_data object containing the
% probability_map data

if isa(mask, 'atlas')
    mask = atlas_get_probability_maps(mask);
end

n_obs2 = size(mask.dat, 2);

% This is done for n images in obj
r = zeros(n_obs2, n_obs);

switch sim_metric
    case 'corr'
        % Correlation
        % Note: There have been some problems with different versions of
        % Matlab having different behaviors for corr.m when they are
        % matrices. The call to corr.m below is replaced with a custom correlation
        % function.
        % r = corr(double(obj.dat), double(mask.dat))';
        
        % vector to matrix correlation formula - anonymous function.
        % a is an N x 1 vector, b is an N x k matrix
        % corr_matrix = @(a, b) ((a-mean(a))' * (b-mean(b)) ./ (length(a) - 1))' ./ (std(a)*std(b)');

        % In practice, the function above would work, but canlab_pattern_similarity
        % excludes missing values (0 or NaN) in data images, which are
        % assumed to be missing. This is a special case for image data, as
        % 0 is assumed to be a missing value.
        % See help canlab_pattern_similarity for more detail.
        
        for im = 1:n_obs2
            
            % r(im, :) = corr_matrix(double(mask.dat(:,im)), double(obj.dat));
            
            r(im, :) = canlab_pattern_similarity(obj.dat, mask.dat(:, im), 'correlation');
            
        end
        
    case 'cosine'
        % Cosine similarity
        
        for im = 1:n_obs2
            %         a = nansum(obj.dat .^ 2) .^ .5; %PK KEEP out of mask for norm
            %         b = nansum(mask.dat(nonemptydat,im ) .^ 2) .^ .5; %PK exlude empty data for norm
            %
            %         r(im, :) = (nansum(bsxfun(@times, obj.dat, mask.dat(:,im))) ./ (a .* b))';
            %
            r(im, :) = canlab_pattern_similarity(obj.dat, mask.dat(:, im), 'cosine_similarity');
        end
        
    case 'overlap'
        % binary overlap
        for im = 1:n_obs2
            r(im, :) = canlab_pattern_similarity(obj.dat, mask.dat(:, im), 'binary_overlap');
        end
        
end

stats.r = r;


if ~doaverage
    
    if ~noplot
        
        switch plotstyle
            case 'wedge'
                % --------------------------------------------------
               
                if ~dofixRange
                    outercircleradius = min(1, max(abs(r)) + .1*max(abs(r)));
                else
                    outercircleradius = fixedrange;
                end
                
                if bicolor
                    % plot negative values in the complementary color
                    %groupColors(2) = {[1 1 1] - groupColors{1}};
                    hh = tor_wedge_plot(r, networknames, 'outer_circle_radius', outercircleradius, 'colors', groupColors, 'nofigure', 'bicolor');
                    
                else
                    
                    hh = tor_wedge_plot(r, networknames, 'outer_circle_radius', outercircleradius, 'colors', groupColors, 'nofigure');
                end
                
            case 'polar'
                % --------------------------------------------------
                
                if ~dofixRange
                    % Plot values for each image in obj
                    [hh, hhfill] = tor_polar_plot({r}, groupColors, {networknames}, 'nonneg');
                else
                    [hh, hhfill] = tor_polar_plot({r}, groupColors, {networknames}, 'nonneg','fixedrange',fixedrange);
                    % Make legend
                    if ~isempty(obj.image_names)
                        han = makelegend(obj.image_names, groupColors);
                    end
                end
                
        end % plotstyle
    end % doplot
    
    % print similarity matrix
    if printTable
        fprintf('\n');
        switch sim_metric
            case 'corr'
                print_matrix(r, {'Name' 'Pearson''s r'}, networknames);
            case 'cosine'
                print_matrix(r, {'Name' 'Cosine Similarity'}, networknames);
            case 'overlap'
                print_matrix(r, {'Name' 'Percent Overlap'}, networknames);
        end
    end
    
% ------------------------------------------------------------
% ------------------------------------------------------------

elseif doaverage
    % Average across replicates (usually participants) and plot means with
    % error regions (shading)
    
    if strcmp(sim_metric,'corr')
        z=fisherz(r'); %transform values
    else
        z=r'; % may need other transformations for other metrics here (SG 2017/09/07)
    end
    
    if exist('compareGroups','var') %if we want to do analysis for multiple groups
        
        groupValues=unique(group);
        g=num2cell(groupValues); %create cell array of group numbers
          
        for i=1:size(z,2) %for each spatial basis do an anova across groups
            
            [p table_group{i} st]=anova1(z(:,i), group, 'off'); %get anova table
            [c,~] = multcompare(st, 'Display', 'off'); %perform multiple comparisons
            multcomp_group{i}=[g(c(:,1)), g(c(:,2)), num2cell(c(:,3:end))]; %format table for output
            
        end
        
        for i=1:size(z,2)
            disp(['Between-group comparisons for ' networknames{i} ':']);
            disp('--------------------------------------');
            disp(['One-way ANOVA: F(' num2str(table_group{i}{2,3}) ','  num2str(table_group{i}{3,3}) ') = ' num2str(table_group{i}{2,5},3) ', P = ' num2str(table_group{i}{2,6},3)])
            disp(' ')
            disp('Multiple comparisons of means:')
            disp(' ');
            print_matrix(cell2mat(multcomp_group{i}), {'Group 1' 'Group 2' 'LCI' 'Estimate' 'UCI' 'P'});
            disp(' ');
        end
          
    else
        group=ones(size(r,2),1); %otherwise all data is from same group
        groupValues=unique(group);
        g=num2cell(groupValues); %creat cell array of group numbers
        
        
    end % Compare groups
    
    
    % Perform test of uniformity for each group
    % ------------------------------------------------------------

    for g=1:length(groupValues)
        
        r_group=r(:,group==groupValues(g));
        z_group=z(group==groupValues(g),:);
        
        stats(g).r = r_group;
        
        % Plot mean and se of values
        m(:,g) = nanmean(r_group')';
        se(:,g) = ste(r_group')';
        
        
        %[h, p, ci, stat] = ttest(r');
        [h, p, ci, stat] = ttest(z_group);
        if strcmp(sim_metric,'corr')
            stats(g).descrip = 'T-test on Fisher''s r to Z transformed point-biserial correlations';
        else
            stats(g).descrip = ['T-test on raw similarity measured by ' sim_metric];  % added SG. See Fisher-Z transformation above.
        end
        stats(g).networknames = networknames;
        stats(g).p = p';
        stats(g).sig = h';
        stats(g).t = stat.tstat';
        stats(g).df = stat.df';
        
        %perform repeated measures anova  (two way anova with subject as the
        %row factor
        [~, stats(g).table_spatial, st]=anova2(z_group(~any(isnan(z_group')),:),1,'off');
        [c,~] = multcompare(st,'Display','off');
        stats(g).multcomp_spatial=[networknames(c(:,1))', networknames(c(:,2))', num2cell(c(:,3:end))];
        
        
        disp(['Table of correlations Group:' num2str(g)]);
        disp('--------------------------------------');
        disp(stats(g).descrip)
        
        print_matrix([m(:,g) stats(g).t stats(g).p stats(g).sig], {'R_avg' 'T' 'P' 'sig'}, networknames);
        disp(' ');
        
    end %groups
    
    % Plot (average + error bars)
    % ------------------------------------------------------------
    if ~noplot
        % groupColors = scn_standard_colors(length(groupValues))'; %
        % removed to enable use of user-defined colors. SG 2017/2/7
        groupColors=repmat(groupColors',3,1);
        groupColors={groupColors{:}};
        
        toplot=[];
        
        for i=1:length(groupValues)
            toplot=[toplot m(:,i)+se(:,i) m(:,i) m(:,i)-se(:,i)];
        end
        
        if dofigure, create_figure('tor_polar'); end
        
        switch plotstyle
            case 'wedge'
                % --------------------------------------------------

                if ~dofixRange
                    outercircleradius = min(1, max(m) + .1*max(m));
                else
                    outercircleradius = fixedrange;
                end
                
                if ~dofixRange
                    outercircleradius = min(1, max(abs(toplot(:))) + .1*max(abs(toplot(:))));
                else
                    outercircleradius = fixedrange;
                end
                
                if bicolor
                    hh = tor_wedge_plot(r_group', networknames, 'outer_circle_radius', outercircleradius, 'colors', groupColors, 'nofigure', 'bicolor');
                else
                    error('bicolor is set to 1 by default - this should not happen. debug me.');
                end
                
            case 'polar'
                % --------------------------------------------------
                
                if ~dofixRange
                    [hh, hhfill] = tor_polar_plot({toplot}, groupColors, {networknames}, 'nonneg');
                else
                    [hh, hhfill] = tor_polar_plot({toplot}, groupColors, {networknames}, 'nonneg', 'fixedrange',fixedrange);
                end
                
                set(hh{1}(1:3:end), 'LineWidth', 1); %'LineStyle', ':', 'LineWidth', 2);
                set(hh{1}(3:3:end), 'LineWidth', 1); %'LineStyle', ':', 'LineWidth', 2);
                
                set(hh{1}(2:3:end), 'LineWidth', 4);
                set(hhfill{1}([3:3:end]), 'FaceAlpha', 1, 'FaceColor', 'w');
                set(hhfill{1}([2:3:end]), 'FaceAlpha', 0);
                set(hhfill{1}([1:3:end]), 'FaceAlpha', .3);
                
                handle_inds=1:3:length(hh{1});
                for g=1:length(groupValues)
                    stats(g).line_handles = hh{1}(handle_inds(g):handle_inds(g)+2);
                    stats(g).fill_handles = hhfill{1}(handle_inds(g):handle_inds(g)+2);
                end
                
                % doaverage
                
                kk = size(toplot, 1);
                mytextsize = 30 ./ (kk.^.3);
                hhtext = findobj(gcf, 'Type', 'text'); set(hhtext, 'FontSize', mytextsize);
                
            otherwise
                error('Unknown plottype');
                
        end % switch plotstyle
    end
    
end % doaverage

% Fill in additional info
for i=1:length(stats)
    stats(i).descrip = 'Rows are networks, columns are input images.';
    stats(i).networknames = networknames;
    stats(i).network_imagenames = imagenames;
    stats(i).inputnames = format_strings_for_legend(obj.image_names);
    stats(i).input_imagenames = obj.image_names;
end
end % function



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Sub-functions
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% % The functions below replaced by load_image_set.  Tor: July 2016

% function [mask, networknames, imagenames] = load_bucknerlab_maps
%
% % Load Bucker Lab 1,000FC masks
% % ------------------------------------------------------------------------
%
% names = load('Bucknerlab_7clusters_SPMAnat_Other_combined_regionnames.mat');
% img = which('rBucknerlab_7clusters_SPMAnat_Other_combined.img');
%
% mask = fmri_data(img, [], 'noverbose');  % loads image with integer coding of networks
%
% networknames = names.rnames(1:7);
% k = length(networknames);
%
% newmaskdat = zeros(size(mask.dat, 1), k);
%
% for i = 1:k  % breaks up into one map per image/network
%
%     wh = mask.dat == i;
%
%     nvox(1, i) = sum(wh);
%
%     newmaskdat(:, i) = double(wh);
%
%
% end
%
% mask.dat = newmaskdat;
%
% imagenames = {img};
% end  % function
%
%
%
%
% function [mask, networknames, imagenames] = load_npsplus
%
% % Load NPS, PINES, Rejection, VPS,
% % ------------------------------------------------------------------------
%
% networknames = {'NPS' 'PINES' 'RomRejPattern' 'VPS'};
%
% imagenames = {'weights_NSF_grouppred_cvpcr.img' ...  % NPS
%     'Rating_Weights_LOSO_2.nii'  ...  % PINES
%     'dpsp_rejection_vs_others_weights_final.nii' ... % rejection
%     'bmrk4_VPS_unthresholded.nii'};
%
% imagenames = check_image_names_get_full_path(imagenames);
%
% mask = fmri_data(imagenames, [], 'noverbose');  % loads images with spatial basis patterns
%
% end  % function
%
%
%
%
%
%
%
% function [mask, networknames, imagenames] = load_kragelemotion
%
% % Load NPS, PINES, Rejection, VPS,
% % ------------------------------------------------------------------------
%
% networknames = {'Amused' 'Angry' 'Content' 'Fearful' 'Neutral' 'Sad' 'Surprised'};
%
% imagenames = { ...
%     'mean_3comp_amused_group_emotion_PLS_beta_BSz_10000it.img' ...
%     'mean_3comp_angry_group_emotion_PLS_beta_BSz_10000it.img' ...
%     'mean_3comp_content_group_emotion_PLS_beta_BSz_10000it.img' ...
%     'mean_3comp_fearful_group_emotion_PLS_beta_BSz_10000it.img' ...
%     'mean_3comp_neutral_group_emotion_PLS_beta_BSz_10000it.img' ...
%     'mean_3comp_sad_group_emotion_PLS_beta_BSz_10000it.img' ...
%     'mean_3comp_surprised_group_emotion_PLS_beta_BSz_10000it.img'};
%
% imagenames = check_image_names_get_full_path(imagenames);
%
% mask = fmri_data(imagenames, [], 'noverbose');  % loads images with spatial basis patterns
%
% end % function
%
%
% function [mask, networknames, imagenames] = load_allengenetics
%
% % Load Allen Brain Atlas project human genetic maps (from Luke Chang)
% % ------------------------------------------------------------------------
%
% networknames = {'5HT' 'Opioid' 'Dopamine' 'NEalpha' 'NEbeta'};
%
% imagenames = { ...
%     'Serotonin.nii' ...
%     'Opioid.nii' ...
%     'Dopamine.nii' ...
%     'AdrenoAlpha.nii' ...
%     'AdrenoBeta.nii' ...
% };
%
% imagenames = check_image_names_get_full_path(imagenames);
%
% mask = fmri_data(imagenames, [], 'noverbose');  % loads images with spatial basis patterns
%
% end % function
%
%
%
%
% function imagenames = check_image_names_get_full_path(imagenames)
%
% for i = 1:length(imagenames)
%
%     if isempty(which(imagenames{i}))
%         fprintf('CANNOT FIND %s \n', imagenames{i})
%         error('Exiting.');
%     end
%
%     imagenames{i} = which(imagenames{i});
% end
% end
%




function [yfit, w, dist_from_hyplane] = crossVal_svm(svmobj,z,y,train_ind,test_ind)

dataobj = data('spider data', z(train_ind,:), y(train_ind,:));


% Training
[res, svmobj] = train(svmobj, dataobj);
res2 = test(svmobj, data('test data', z(test_ind,:), []));
yfit = res2.X;



if size(y,2)==1
    w = get_w(svmobj)';
    b0 = svmobj.b0;
    dist_from_hyplane = z(test_ind,:) * w + b0;
else
    for i = 1:size(res2.X,2)
        w(:,i) = get_w(svmobj{i})';
        b0(i) = svmobj{i}.b0;
        dist_from_hyplane(:,i) = z(test_ind,:) * w(:,i) + b0(i);
        
    end
end
end