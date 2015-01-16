% function [classes, stats] = cluster_princomp2(y, outcome, covariates, varargin)
%
% tor wager, in progress April 2007
%
% Example: Predict a behavioral outcome from data in a set of clusters,
% with covariates
% -----------------------------------------------------------------------
% for i = 1:length(cl), y(:,i) = cl(i).CONTRAST.indiv_data; end
% covariates = EXPT.cov(:,2:end);
% outcome = EXPT.cov(:,1);
% covnames = EXPT.covnames(2:end);
% [classes, stats] = cluster_princomp2(y, outcome, covariates, covnames, 'cl', cl);
%
% With mediation SETUP file, on rank data
% This way removes all covariates first, so they cannot be related to the
% outcome in the end.
% c = []; c.covs_nointerest = SETUP.X; c.outcome = SETUP.Y; doranks = 1;
% c = nmdsfig_tools('get_data_matrix',c,cl,'timeseries',1,[],doranks);
% [classes, stats] = cluster_princomp2(c.dat, c.outcome, c.covs_nointerest, covnames, 'cl', cl);

function [classes, stats] = cluster_princomp2(y, outcome, covariates, covnames, varargin)
    
% defaults
% -----------------------------------------------------------------------
niter = 1000;

alph = .05;     % for stepwise regression
doplot = 1;
dosave = 0;

%cnames = {'red' 'blue' 'yellow' 'green' 'orange' 'cyan' 'purple'};
%colors = {[1 0 0] [0 0 1] [1 1 0] [0 1 0] [1 .5 0] [0 .5 1] [.5 0 1]};
cnames = {'yellow' 'green' 'blue' 'red' 'cyan' 'magenta' 'black'};
colors = {'yo' 'go' 'bo' 'ro' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};



for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'cl', cl = varargin{i+1};
            case 'plot', doplot = varargin{i+1};
            case 'save', dosave = 1;

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


nanvals = any(isnan(y));
if any(nanvals)
    fprintf(1,'Removing columns with NaNs : %3.0f \n', find(nanvals));
    
    y(:, nanvals) = [];
    if exist('cl', 'var'), cl(nanvals) = []; end
    
end

[nobs, nregions] = size(y);

% principal components
% -----------------------------------------------------------------------

[pc, stats] = pca_npm(y, niter);
ncomps = sum(stats.sig);

scn_export_papersetup(400);
set(gcf,'Tag', 'eigenvalues')

if dosave
    saveas(gcf,'prediction_pca','png');
end
% scores = stats.score(:, stats.sig);

% classification into classes
% -----------------------------------------------------------------------
maxtotest = min(7, nregions);
[bestpval,classes]=testclustnew(pc,2:maxtotest,ncomps,niter,[],'keep','average');
nclasses = max(classes);

scn_export_papersetup(400);
set(gcf,'Tag', 'clustering')

if dosave
    saveas(gcf,'prediction_clustering','png');
end

% % for i = 2:nregions
% %     classes = clusterdata(pc, 'linkage', 'average', 'maxclust', i);
% %     s = silhouette(pc, classes); sil(i) = mean(s);
% % end
%nclasses = size(stats.class, 2);


class_averages = zeros(nobs, nclasses);
for i = 1:nclasses
    wh = (classes == i);
    class_averages(:, i) = mean(y(:, wh), 2);
end

%STEP = stepwise_tor(stats.classdata(:,1:3), EXPT.cov(:,1));

% names 
% -----------------------------------------------------------------------
names = cell(1, nclasses);
for i = 1:nclasses
    names{i} = ['Brain ' num2str(i) ' (' cnames{i} ')'];
end

names = [names covnames];

% model: predictors
% -----------------------------------------------------------------------
X = [class_averages covariates];

% stepwise regression
% -----------------------------------------------------------------------

STEP = stepwise_tor(X, outcome, names, alph);

stats.classes = classes;
stats.class_averages = class_averages;
stats.names = names;

wh = logical(STEP.inmodel);
X = X(:, wh);
names = names(:, wh);
stats.STEP = STEP;
stats.sig_X = X;
stats.sig_names = names;
stats.outcome = outcome;
stats.covariates = covariates;

% plot
% -----------------------------------------------------------------------
if doplot && exist('cl','var')
    domontage = 0;
    cluster_orthviews_classes(cl,classes,[],'axial',domontage);
    spm_orthviews('Xhairs', 'on')
    
    scn_export_papersetup(600);
    set(gcf,'Tag', 'slices')
    if dosave
        saveas(gcf,'prediction_slices','png');
    end
    
end

nsig = size(X,2);

if doplot && nsig > 0
    
    
    tor_fig(1, nsig); 
    set(gcf,'Tag','scatterplots');
    scn_export_papersetup(300);
    
    for i = 1:nsig
        subplot(1, nsig, i);
    
        prplot(stats.outcome,stats.sig_X,i,0,colors(i));
        title(names{i})
        xlabel(' ');
        ylabel(' ');
        
    end
    
    if dosave
        saveas(gcf,'prediction_scatterplots','png');
        
        save cluster_princomp2_stats stats
        
        print_output;
    end
    
end
    
end

% see cluster_orthviews_classes
% % 
% % addflag = 'newfig';
% % 
% % for i = 1:nclasses
% %     
% %     % regions in this group
% %     wh = classes == i;  %logical(stats.class(:,i)');
% %     
% %     cluster_orthviews(cl(wh), colors(i), 'overlay',[], addflag, 'solid');
% %     addflag = 'add';
% %     
% %     
% % end

function print_output

diary cluster_princomp2_output.txt

fprintf('\n\ncluster_princomp2\n------------------------------\n\n');
fprintf('\nObservations: %3.0f', nobs);
fprintf('\nRegions (var): %3.0f', nregions);
fprintf('\nIterations: %3.0f', niter);
fprintf('\nSignificant PCs: %3.0f', ncomps);
fprintf('\nSignificant classes: %3.0f', nclasses);
fprintf('\nSig. stepwise predictors: %3.0f\n', nsig);

fprintf('\n')
STEP = stepwise_tor(X, outcome, names, alph);
STEP;

if exist('cl', 'var') && isfield(cl, 'shorttitle')
    for i = 1:nclasses
        fprintf('\nClass %3.0f\n', i);
        for j = find(classes == i)
            fprintf('%s  ', cl(j).shorttitle);
        end
        fprintf('\n')
    end
end

fprintf('\n------------------------------\n');

diary off

end


