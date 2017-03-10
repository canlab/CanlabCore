function [c,cl] = cluster_nmdsfig(cl,fldname,colm,varargin)
% function [c,cl] = cluster_nmdsfig(cl,fldname,colm,[behavior],[ancova var],[covs of no interest],[overlay],[niter])
%
% Makes a multidim scaling figure from clusters
% Plus some other plots
%
% Example:
% c = cluster_nmdsfig(cl,'BARPLOT.data',2);
%
% c is a structure with lots of info
% cl has names appended, if not entered previously
%
% behavioral input is of three types:
% outcome       covariate of interest (1 only)
% ancova_var        categorical ANCOVA coding var (2-levels, one only);
%                   generally not of interest
%                   IF ancova_var is entered, "PPI" of group x brain
%                   covariance interaction will be computed
% covs of no int    of no interest; will be removed from brain and behavior
%                   before performing ANCOVA (for simplicity); 1 or more OK
%
% Tor Wager
%
% Examples:
%
% Use average data from cl.timeseries and correlate with EXPT.behavior.
% Covary out effects of order before running models; no ANCOVA
% [c,cl] = cluster_nmdsfig(cl,'timeseries',1,EXPT.behavior',[],EXPT.order2);
%
% More examples:
% c =
% cluster_nmdsfig(cl,'CONTRAST.yadj',1,repmat(EXPT.cov(:,1),4,1),[],[],EXPT.overlay,1000);
% [c,cl] = cluster_nmdsfig(cl,'CONTRAST.data',1,R.X(:,1),[],R.X(:,2),EXPT.overlay,200);
% [c,cl] = cluster_nmdsfig(cl,'timeseries',1,SETUP.Y,[],SETUP.X,EXPT.overlay,1000);

% subdirectory: prompt
dosubdir = input('Create a subdirectory for png images? (type name or return to use current dir) ','s');
if ~isempty(dosubdir)
    mkdir(dosubdir)
    cd(dosubdir)
end

diary cluster_nmdsfig_output.txt
doranks = 1;
addbtwn = 0;

if length(varargin) > 0, c.outcome = varargin{1}; else  c.outcome = [];end
if length(varargin) > 1, c.ancova_codes = varargin{2}; else  c.ancova_codes = [];end
if length(varargin) > 2, c.covs_nointerest = varargin{3}; else  c.covs_nointerest = [];end
if length(varargin) > 3, overlay = varargin{4}; else  overlay = [];end
if length(varargin) > 4, niter = varargin{5}; else  niter = 1000;end

cl(1).outcome = c.outcome; cl(1).ancova_codes = c.ancova_codes; cl(1).covs_nointerest=c.covs_nointerest;

c.X = [c.outcome, c.covs_nointerest];
contrasts = [];

c.colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};

if ~isfield(c,'btwn_names') || isempty(c.btwn_names)
    c.btwn_names = {'Behav'};
end


% -----------------------------------------------------------------------
% get data from selected field
% remove covariates of no interest from predictors and data
% -----------------------------------------------------------------------

c = nmdsfig_tools('get_data_matrix',c,cl,fldname,colm,contrasts,doranks);

if addbtwn
    % add between-subjects behavioral scores to dat; just to get stim coords; remove later...
    c.dat = [c.outcome c.dat];
end

% -----------------------------------------------------------------------
% get names
% -----------------------------------------------------------------------

cl = cluster_names(cl,1);
for i = 1:length(cl), c.names{i} = cl(i).shorttitle; end


% -----------------------------------------------------------------------
% get correlations, significance, and number of dimensions
% If groups are entered and ANCOVA is specified, handle that
% Also get c.D, dissimilarity matrix for MDS
% -----------------------------------------------------------------------

c = nmdsfig_tools('get_correlations',c);


% MDS decomposition and plotting
% -----------------------------------------------------------------------
% get number of dimensions that carry most of the variance.
% This is to simplify the space and make it less sparse --
% expected to give better clustering results.
%
% Reduce to distances in ndims space, check this against actual distances
% -----------------------------------------------------------------------

% graphic check on number of dimensions, and MDS
[c.GroupSpace,c.obs,c.implied_dissim] = shepardplot(c.D,[]);
drawnow
c.ndims = size(c.GroupSpace,2); % set here; graphical decision within shepardplot

nclust = 2:(size(c.GroupSpace,1) / 2)+1;     % choose number of clusters to test

save_figure('cluster_nmdsfig_shepardplot');

if addbtwn
    % get coordinates for behavior along with the others
    fprintf(1,'Using between-subjects covariates in MDS model.\n');

    % save coordinates in MDS space for btwn data
    c.btwn_coords = c.GroupSpace(1,:);

    % remove from GroupSpace
    c.GroupSpace = c.GroupSpace(2:end,:);

    % remove from r
    c.r = c.r(2:end,2:end);
end

% -----------------------------------------------------------------------
% testcluster:
%
% Permute objects in space, get null hypothesis cluster quality, and test
% observed cluster solution against this.
% -----------------------------------------------------------------------

% pval is p-values for each number of clusters tested
% classes = class ("network") assignments for best clustering solution
[c.ClusterSolution.pvals,c.ClusterSolution.classes, c.ClusterSolution.classnames ...
    c.ClusterSolution.X,c.ClusterSolution.include,c.ClusterSolution.names]= ...
    testclustnew(c.GroupSpace,nclust,c.ndims,niter,c.names,'keep','average');

if ~any(c.ClusterSolution.pvals < .05)
    % nothing significant
    c.ClusterSolution.classes = ones(size(c.ClusterSolution.classes));
end

%[bestpval,bestXc,bestnames,bestX,where,clustnames]=testcluster(X,clust,[r],[nperm],[names],[remove],[linkagetype]);

save_figure('cluster_nmdsfig_testcluster');



% -----------------------------------------------------------------------
% factor scores with reduced dimensions
% -----------------------------------------------------------------------
r = c.r;
[tmp c.eigenvalues] = cmdscale(c.D);
sqL = sqrt(diag(c.eigenvalues(1:c.ndims)));
B = inv(r) * c.GroupSpace(:,1:c.ndims) * sqL;   % inv(r) * factor loading matrix A = VsqrtL
c.factor_scores = c.dat * B;

% -----------------------------------------------------------------------
% Stepwise regression of dimensions on behavior
% -----------------------------------------------------------------------

if ~isempty(c.outcome)
    fprintf(1,'\n');
    fprintf(1,'Prediction of behavior with component scores (from classic MDS)\n');

    DATA.INDSCAL.STEPWISE = stepwise_tor(c.factor_scores,c.outcome);
    fprintf(1,'\n');
end



% -----------------------------------------------------------------------
% % plot nmds figures of this
% -----------------------------------------------------------------------
dosideplots = 0;
figname = nmdsfig_tools('nmdsfig_plot',c, addbtwn, dosideplots, 'nofill');

save_figure(figname);


% -----------------------------------------------------------------------
% separate plots of high/low behavior
% -----------------------------------------------------------------------
if ~isempty(c.outcome)  % plot interactions between beh and regional covariance
    create_wide_figure;
    subplot(1,2,1);
    [r,p]=corrcoef(c.dat(c.outcome>median(c.outcome),:)); r(r>.9999) = 1; sig = sign(r) .* (p < .05);
    nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',sig);
    title('High behavior')

    subplot(1,2,2)
    [r,p]=corrcoef(c.dat(c.outcome<=median(c.outcome),:)); r(r>.9999) = 1; sig = sign(r) .* (p < .05);
    nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',sig);
    title('Low behavior')
    drawnow
end

if ~isempty(c.ancova_codes)  % plot interactions between beh and regional covariance
    tor_fig(1,3); subplot(1,3,1);
    [r,p]=corrcoef(c.dat(c.ancova_codes>median(c.ancova_codes),:)); r(r>.9999) = 1; sig = sign(r) .* (p < .05);
    nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',sig);
    title('High ANCOVA group')

    subplot(1,3,2)
    [r,p]=corrcoef(c.dat(c.ancova_codes<median(c.ancova_codes),:)); r(r>.9999) = 1; sig = sign(r) .* (p < .05);
    nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',sig);
    title('Low ANCOVA group')
    drawnow

    subplot(1,3,3)
    [r,p]=corrcoef(c.dat(c.ancova_codes<median(c.ancova_codes),:)); r(r>.9999) = 1; sig = sign(r) .* (p < .05);
    nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',c.STATS.siginteract);
    title('Interactions')


end


% -----------------------------------------------------------------------
% Optional plots of interactions and correlations
% -----------------------------------------------------------------------
if c.doancova && ~isempty(c.ancova_codes)
    [i,j] = find(triu(c.STATS.siginteract));
    if ~isempty(i),
        go = input([num2str(length(i)) ' significant behavior x covariance interactions.  Plot? ']);
        if go
            for ind = 1:length(i)
                cluster_orthviews(cl(i(ind)),{[1 0 0]});
                cluster_orthviews(cl(j(ind)),{[0 0 1]},'add');
                h = axes('Position',[.55 .05 .42 .42]);
                [b,t,p] = ancova(c.ancova_codes,c.dat(:,i(ind)),c.dat(:,j(ind)),1);
                xlabel(cl(i(ind)).shorttitle); ylabel(cl(j(ind)).shorttitle);
                title([]);
                input('Save and press a key for the next one')
            end
        end
    end
end

if c.doancova, c.doancova = 'Yes'; else c.doancova = 'No'; end

if ~any(c.ClusterSolution.pvals < .05)
    go  = input('No significant clustering.  Proceed with analyses anyway?');
    if ~go, return, end
end





% -----------------------------------------------------------------------
% Get average correlations within and between clusters
%
% Apply this clustering solution to the clusters
% Get individual 'network' scores and use them as predictors of behavior
% in multiple regression
% -----------------------------------------------------------------------

c = nmdsfig_tools('apply_clusters',c);

c = nmdsfig_tools('predict_behavior',c,'regions');
c = nmdsfig_tools('predict_behavior',c,'classes');




% -----------------------------------------------------------------------
% brain plots
% -----------------------------------------------------------------------
if isempty(overlay)
    cont = input('Show brain slices? ');
    if ~cont
        diary off

        save nmds_output c cl
        cd ..
        return
    end

    overlay = spm_get(1,'*img','Select overlay image');
else
    % we have overlay, implies we want slices

    cluster_orthviews_classes(cl,c.ClusterSolution.classes,overlay,'axial',0);
    save_figure('cluster_slices_axial');

    %     cluster_orthviews_showcenters(cl,'coronal',overlay);
    %     save_figure('cluster_slices_coronal');
    %
    %     cluster_orthviews_showcenters(cl,'saggital',overlay);
    %     save_figure('cluster_slices_saggital');
end

diary off

save nmds_output c cl
cd ..

return




% -----------------------------------------------------------------------

% -----------------------------------------------------------------------


% Sub-functions


% -----------------------------------------------------------------------

% -----------------------------------------------------------------------






% -----------------------------------------------------------------------
% Utility functions
% -----------------------------------------------------------------------
function save_figure(myname)
f1 = gcf;
scn_export_papersetup;
saveas(f1,myname,'png');
return


function f1 = create_wide_figure
scnsize = get(0,'ScreenSize');
f1 = figure('position',[50 50 scnsize(3)-100 scnsize(4)/2],'color','white');
return



function [pthr,sig] = fdr_correct_pvals(p,r)

psq = p; psq(find(eye(size(p,1)))) = 0;
psq = squareform(psq);
pthr = FDR(p,.05);
if isempty(pthr), pthr = 0; end

sig = sign(r) .* (p < pthr);

return


