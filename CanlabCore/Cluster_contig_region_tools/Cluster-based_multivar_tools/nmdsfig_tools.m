%
% Tools for using nmds multidimensional scaling plots and data structure
% the data structure variable in this code is called c (arbitrarily)
% See cluster_nmdsfig for a function that uses many of these.
%
% By Tor Wager, Aug-Sept 2006
% V 1.2, April 2007
%
% ---------------------------------------------------------------------
% --
% Get data matrix, rank [optional], adjust for nuisance [optional]
% -----------------------------------------------------------------------
% c = nmdsfig_tools('get_data_matrix',c,cl,fldname,colm,contrasts,doranks)
% c = nmdsfig_tools('get_data_matrix',c,cl,'CONTRAST.avg',1,[],1);
%
% define input c.covs_nointerest with covariates to remove them.
% define input c.outcome to also remove covs from outcome.
%
% -----------------------------------------------------------------------
% get correlations, significance, and number of dimensions
% If groups are entered and ANCOVA is specified, handle that
% -----------------------------------------------------------------------
% c = nmdsfig_tools('get_correlations',c);
%
% if c.ancova_codes is input, does ancova analysis and returns stats
%
% -----------------------------------------------------------------------
% get cluster solution for a set of data (objectsByDimsData or
% GroupSpace)
% -----------------------------------------------------------------------
% c = cluster_solution(c, objectsByDimsData, n_cluster_range, niter, names)
%
% names is cell of one name per object
% to set up input, use get_correlations, and then this:
% [c.GroupSpace,c.obs,c.implied_dissim] = shepardplot(c.D,[]);
% c = nmdsfig_tools('cluster_solution',c, c.GroupSpace, 2:5, 1000, []);
%
% -----------------------------------------------------------------------
% Make a connectivity figure; with optional add btwn vars and side plots
% Uses stats from correlations
% -----------------------------------------------------------------------
% nmdsfig_tools('nmdsfig_plot',c, addbtwn, dosideplots, fillstr)
% addbtwn, dosideplots are 1 or 0, fillstr is 'fill' or something else
% uses colors in c.colors if it exists.
% -----------------------------------------------------------------------
% Make a connectivity figure; WITH outcome var and covariates on figure
% -----------------------------------------------------------------------
% DO THIS USING c STRUCTURE AFTER CLUSTER SOLUTION IS DONE
% nmdsfig_tools('nmdsfig_plot_withcovs',c, cl, dorank)
% needs fields: c.outcome c.covs_nointerest
%
% -----------------------------------------------------------------------
% Apply a classification to a set of data and r matrix, get class avgs
% -----------------------------------------------------------------------
% c = nmdsfig_tools('apply_clusters',c)
% c.class_avg_dat contains class averages
% needs c.ClusterSolution.classes, c.r
%
% Can turn off interactive naming
% c = nmdsfig_tools('apply_clusters',c, 0);
%
% -----------------------------------------------------------------------
% Get class avgs from some other data matrix using an existing c struct
% -----------------------------------------------------------------------
% class_avg_dat = nmdsfig_tools('get_class_avgs',c,dat)
% needs c.ClusterSolution.classes, dat (data matrix)
%
% -----------------------------------------------------------------------
% Predict behavior from either class average data or regions
% -----------------------------------------------------------------------
% c = nmdsfig_tools('predict_behavior',c,'classes');
% c = nmdsfig_tools('predict_behavior',c,'regions');
%
% -----------------------------------------------------------------------
% Add lines to a graph
% -----------------------------------------------------------------------
% handles = nmdsfig_tools('drawlines',coords,significance_mtx, bendval);
% specify pos/neg colors and styles:
% nmdsfig_tools('drawlines',c.GroupSpace,c_compare.sig,[1 0 0;0 0 1],{':',':'}, .05);
%
% -----------------------------------------------------------------------
% % Remove lines from a graph
% -----------------------------------------------------------------------
% nmdsfig_tools('removelines');
%
% -----------------------------------------------------------------------
% Add stars where stepwise params are significant
% -----------------------------------------------------------------------
% nmdsfig_tools('stepwise_stars',c,'classes')
% nmdsfig_tools('stepwise_stars',c,'regions')
%
% or where robust regression with behavior is significant
% nmdsfig_tools('robust_correl_stars',c,wh_reg)
%
% -----------------------------------------------------------------------
% Get unique colors for each object, but similar colors within a class
% -----------------------------------------------------------------------
% colors = nmdsfig_tools('class_colors',classes,basecolors)
%
%
% -----------------------------------------------------------------------
% Connect two points x and y in 3-D space with a straight or curved line
% -----------------------------------------------------------------------
% out = nmdsfig_tools('connect3d',x,y, varargin)
% x and y are [x,y,z] coordinate triples, or clusters/region objects with mm_center field.
% out.h is handles, out.coords are also returned
%
% Optional inputs: Keyword followed by value:
% color,thickness,bendpercent, nstreamlines, streamlineshift, nsamples
% bendpercent can be scalar or [x y z] triple for bend in each direction
% streamlineshift is multiplier for variation in x y z loc across streamlines
%
% Note: for 3-D graphs, use cluster_nmdsfig_glassbrain.m
%
% -----------------------------------------------------------------------
% Connect two points x and y in 2-D space with a straight or curved line
% -----------------------------------------------------------------------
% out = nmdsfig_tools('connect2d',x,y,color,thickness,bendpercent)
% x and y are [x,y] coordinate doubles, or clusters with mm_center field.
% out.h is handles, out.coords are also returned
% bendpercent can be scalar or [x y] pair for bend in each direction
% e.g.,
% out = nmdsfig_tools('connect2d', [-.3 0], [-.1 .1], 'b', 3, .1);
%
% SEQUENCES FOR ANALYSIS
% and examples
% -----------------------------------------------------------------------
% See also: cluster_nmdsfig, for a batch sequence
% -----------------------------------------------------------------------
% Use average data across conditions to get basic structure:
% c = nmdsfig_tools('get_data_matrix',c,cl,'CONTRAST.avg',1,[],1);
% c = nmdsfig_tools('get_correlations',c);
% nmdsfig_tools('nmdsfig_plot',c, 0, 0);
% c = nmdsfig_tools('apply_clusters',c);
% c = nmdsfig_tools('predict_behavior',c,'regions');
% c = nmdsfig_tools('predict_behavior',c,'classes');
%
% Get new data and impose the class structure on it:
% c = nmdsfig_tools('get_data_matrix',c,cl,'CONTRAST.data',1:2,[1 -1],0);
% c = nmdsfig_tools('predict_behavior',c,'regions');
% c = nmdsfig_tools('get_correlations',c);
% c = nmdsfig_tools('apply_clusters',c);
% c = nmdsfig_tools('predict_behavior',c,'classes')
%
% Create a basic figure based on the averages,
% and then use new data to get lines
% add stars where there are sig. correls with behavior
% c = nmdsfig_tools('get_data_matrix',c,cl,'CONTRAST.avg',1,[],1);  % ranks, for robustness
% c = nmdsfig_tools('get_correlations',c);
%f1 = nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'sig',c.STATS.sigmat,'names',c.names,'sizescale',[4 16]);
% nmdsfig_tools('removelines');
% c = nmdsfig_tools('get_data_matrix',c,cl,'CONTRAST.data',1:2,[1 -1],0);
% c = nmdsfig_tools('get_correlations',c);
% nmdsfig_tools('drawlines',c.GroupSpace,c.STATS.sigmat);
% nmdsfig_tools('robust_correl_stars',c,1);
%
% Get correl diffs among class averages
% get_data_matrix adjusts for covs of no interest, ranks, and centers
%c = nmdsfig_tools('get_data_matrix',c,cl,'CONTRAST.data',2,[],0); yON = c.dat; % get rank data
%c = nmdsfig_tools('get_data_matrix',c,cl,'CONTRAST.data',1,[],0); yOFF = c.dat; % get rank data
%yONavgs = nmdsfig_tools('get_class_avgs',c,yON);
%yOFFavgs = nmdsfig_tools('get_class_avgs',c,yOFF);
% ycomp_struct = correl_compare_dep(yONavgs,yOFFavgs,.05,0); % data already
% %ranked, so no need to redo
%
%
% Complete example for pre-appraisal mediation data:
% For set-up, see mediation_brain, mediation_brain_results and
% mediation_brain_results_detail
%
%     c = []; c.covs_nointerest = [SETUP.X SETUP.covariates]; c.outcome = SETUP.Y; doranks = 1;
%     for i = 1:length(cl), c.names{i} = cl(i).shorttitle; end
%     c = nmdsfig_tools('get_data_matrix',c,cl,'timeseries',1,[],doranks);
%     c = nmdsfig_tools('get_correlations',c);
%     [c.GroupSpace,c.obs,c.implied_dissim] = shepardplot(c.D,[]);
%     c = nmdsfig_tools('cluster_solution',c, c.GroupSpace, 2:5, 1000, []);
%     c = nmdsfig_tools('apply_clusters',c);
%     c = nmdsfig_tools('predict_behavior',c,'classes');
%     cluster_orthviews_classes(cl,c.ClusterSolution.classes, EXPT.overlay, 'sagittal', 0);
%     nmdsfig_tools('nmdsfig_plot',c, 0, 0, 'fill');
% c_complete = nmdsfig_tools('nmdsfig_plot_withcovs',c, cl, doranks)

function [c, varargout] = nmdsfig_tools(meth,varargin)

c = [];

% First define required inputs
% -----------------------------------------------------------------------
switch meth
    
    case 'get_data_matrix'
        % This should execute the function:
        % c = get_data_matrix(c,cl,fldname,colm,contrasts,doranks)
        innames = {'c' 'cl' 'fldname' 'colm' 'contrasts' 'doranks'};
        
    case 'nmdsfig_plot'
        % nmdsfig_plot(c, addbtwn, dosideplots)
        innames = {'c' 'addbtwn' 'dosideplots' 'fillstr'};
        
    case 'nmdsfig_plot_withcovs'
        innames = {'c' 'cl' 'doranks'};
        
        
    case {'apply_clusters','get_correlations'}
        % c = apply_clusters(c)
        innames = {'c', 'dointeractive', 'corrtype'};
        if length(varargin) < 2, varargin{2} = []; end % optional; do interactive 1/0
        if length(varargin) < 3, varargin{3} = 'r'; end % optional; corrtype, 'r' = Pearson't (default)
        
    case {'cluster_solution'}
        innames = {'c' 'objectsByDimsData' 'n_cluster_range' 'niter' 'names'};
        
        % c = cluster_solution(c, objectsByDimsData, n_cluster_range, niter, names);
        
    case 'get_class_avgs'
        % class_avg_dat = nmdsfig_tools('get_class_avgs',c,dat)
        innames = {'c' 'dat'};
        
    case 'predict_behavior'
        %c = predict_behavior(c,meth)
        innames = {'c' 'predmeth'};
        
    case 'drawlines'
        % [hhp,hhn] = drawlines(pc,sigmat,sigcol,legmat,[color],[style])
        
        innames = {'pc','sigmat','color','style','bendval'};
        % optional: color/style format strings
        if length(varargin) < 3, varargin{3} = []; end % optional
        if length(varargin) < 4, varargin{4} = []; end % optional
        if length(varargin) < 5, varargin{5} = []; end % optional
        if length(varargin) > 5, error('Too many input arguments.'); end
        
    case 'removelines'
        innames = {};
        
    case 'stepwise_stars'
        innames = {'c' 'meth'};
        
    case 'robust_correl_stars'
        innames = {'c' 'wh_reg'};
        
    case 'class_colors'
        innames = {'classes' 'basecolors'};
        if length(varargin) < 2, varargin{2} = []; end % optional
        
    case 'connect3d'
        % h = nmdsfig_tools('connect3d',x,y,color,thickness,bendpercent,[nsamples])
        % innames = {'x' 'y' 'color','thickness','bendpercent','nsamples'};
        % This has been updated to use inputParser object (2019)
        innames = varargin;
        %if length(varargin) < 6, varargin{6} = []; end % optional
        % Skip all the stuff below
        c = connect3d(varargin{:});
        return
        
    case 'connect2d'
        % out = connect2d(x,y,color,thickness,bendpercent,varargin)
        innames = {'x' 'y' 'color','thickness','bendpercent','nsamples'};
        if length(varargin) < 6, varargin{6} = []; end % optional
        
    otherwise
        error('Unknown method string.')
end

% this generic code block reads varargin into variables for
% do it this way partly to keep track of vars by name
% -----------------------------------------------------------------------
if length(varargin) < length(innames)
    error('Wrong number of input arguments.  Check help on method you are using.')
end
for i = 1:length(innames)
    eval([innames{i} ' = varargin{' num2str(i) '};']);
end

% Now execute subfunction
% -----------------------------------------------------------------------
str = ['c = ' meth '('];
for i = 1:length(varargin)
    if i ~= 1, str = [str ',']; end
    str = [str innames{i}];
end
str = [str ');'];

%disp(str)
eval(str)


return









% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
%
%
% Sub-functions
%
%
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------




% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% get data from selected field
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

function c = get_data_matrix(c,cl,fldname,colm,contrasts,doranks)

if ~isfield(c, 'outcome')
    c.outcome = [];
end

if ~isfield(c, 'covs_nointerest')
    c.covs_nointerest = [];
end

dat = [];
for i = 1:length(cl)
    str = ['dat1 = cl(i).' fldname '(:,[' num2str(colm) ']);'];
    eval(str)
    
    if ~isempty(contrasts), dat1 = dat1 * contrasts'; end
    
    if doranks, dat1 = rankdata(dat1); end
    dat=[dat dat1];
end

if doranks
    disp('Returning rank data in c.dat');
    c.dat_descrip = 'Dat and outcome data: Ranked';
    
    if ~isempty(c.outcome)
        disp('Ranking outcome data.')
        c.outcome = rankdata(c.outcome);
    end
    
    if ~isempty(c.covs_nointerest)
        disp('Ranking covariate data.')
        for j = 1:size(c.covs_nointerest, 2)
            c.covs_nointerest(:,j) = rankdata(c.covs_nointerest(:,j));
        end
    end
    
    
else
    disp('Returning data in c.dat');
    c.dat_descrip = 'Dat and outcome data: Original (not ranked) ';
end

c.doranks = doranks;
%if length(varargin) > 0,dat = [c.outcome dat];end

% remove covariates of no interest, if any
% if data is not rank data, use robust regression
% otherwise, use OLS
% -----------------------------------------------------------------------

if isempty(c.covs_nointerest)
    % do nothing
    disp('nmdsfig_tools get_data_matrix : no covariates entered.')
    c.dat_descrip = [c.dat_descrip ' raw (not adjusted)'];
else
    disp('nmdsfig_tools get_data_matrix : Removing covariate(s) of no interest (using IRLS).')
    c.dat_descrip = [c.dat_descrip ' Adjusted for covariates with IRLS'];
    if ~isfield(c, 'outcome'), c.outcome = []; end
    
    if ~isempty(c.outcome)
        % partialcorr gets adjusted x, y, and correls; here just use it to
        % get adjusted x
        %c.outcome  = partialcor([c.outcome c.covs_nointerest],ones(size(c.outcome,1),1),1);
        c.outcome =  adjust_y(c.covs_nointerest,c.outcome,doranks);
    end
    
    for i = 1:size(dat,2)
        %dat(:,i) = partialcor([dat(:,i) c.covs_nointerest],ones(size(dat,1),1),1);
        [dat(:,i), b(:, i)] =  adjust_y(c.covs_nointerest,dat(:,i),doranks);
    end
    
    c.nuisance_betas = b;
end

c.dat = dat;


return

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% % Adjust x and y for regressors of no interest

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
function [r, b] = adjust_y(x,y,doranks)

if doranks
    X = [x ones(size(x,1),1)];
    b = pinv(X) * y;
    r = y - X * b;
else
    [b,stats]=robustfit(x,y);
    r = stats.resid;
end

return





% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Make an NMDSfig figure from a c structure

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

function figname = nmdsfig_plot(c, addbtwn, dosideplots, fillstr)

% Required fields: (See cluster_nmdsfig, for example)
% requires c. ...
% STATS.sigmat/sigmat2/siginteract
% outcome
% ClusterSolution.classes
% names
% GroupSpace (stim. coords)
% X
% r
%
% btwn_coords (for addbtwn option)
% factor_scores (for sideplots option)

f1 = create_figure;

colors = {'yo'  'bo' 'go' 'ro' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};

if isfield(c,'colors'), colors = c.colors; end


% for side plots
%subplot('position',[0.3 0.3 0.5 0.6]);
if ~isfield(c,'outcome'), c.outcome = []; end

if ~isempty(c.outcome)  % plot interactions between beh and regional covariance
    if isfield(c.STATS,'siginteract')
        % we have a group x cov interaction in ANCOVA
        nmdsfig(GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',c.STATS.sigmat, ...
            'legend',{'Pos' 'Neg' 'Pos x Beh' 'Neg x Beh'},'sig2',c.STATS.siginteract, ...
            'colors', colors, fillstr);
        figname = 'cluster_nmdsfig_ancova_interact';
    else
        nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',c.STATS.sigmat,'sig2', ...
            c.STATS.sigmat2,'legend',{'Pos' 'Neg'},'sizescale',[4 16], ...
            'colors', colors, fillstr);
        figname = 'cluster_nmdsfig_mdsfig';
        
    end
else
    nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',c.STATS.sigmat,'sig2', ...
        c.STATS.sigmat2,'legend',{'Pos' 'Neg'},'sizescale',[4 16], ...
        'colors', colors, fillstr);
    figname = 'cluster_nmdsfig_mdsfig';
end

try
    nmdsfig_legend(c.GroupSpace,c.r);
catch
end

save_figure('cluster_nmdsfig_legend',550);

% behavior
if addbtwn
    if ~isfield(c, 'btwn_coords'), error('You must enter c.btwn_coords to use btwn-subjects breakdown plot.  See cluster_nmdsfig.'); end
    
    hold on; plot(c.btwn_coords(1),c.btwn_coords(2),'ko','MarkerSize',12,'LineWidth',3);
end


% MDSfig sideplots
% -----------------------------------------------------------------------

dosideplots = 0;

if dosideplots
    
    % X-axis subplot: Component 1
    subplot('position',[0.25 0.1 0.6 0.1]);     hold on;
    set(gca,'FontSize',16);
    fact = c.factor_scores(:,1);
    if ~isempty(c.outcome),
        fact=[fact c.outcome]; fact=sortrows(fact,2);
        xlab='Activation score for each participant: Component 1';
        
        lo = find(fact(:,2)<median(fact(:,2)));
        hi = find(fact(:,2)>median(fact(:,2)));
        
        plot(lo,fact(lo,1),'ko','LineWidth',2);
        plot(hi,fact(hi,1),'k^','LineWidth',2);
        set(gca,'XTick',[1 max(lo)+.5 size(fact,1)],'XTickLabel',{'Lowest' 'Median' 'Highest'})
        title(xlab)
        xlabel('Behavior');
        ylabel('Score')
    else
        xlab='Scores on Component 1';
        plot(fact(:,1),'ko','LineWidth',2);
        title(xlab); xlabel('Index'); ylabel('Score');
    end
    
    % Y-axis subplot: Component 2
    subplot('position',[0.1 0.3 0.1 0.6]);   hold on;
    set(gca,'FontSize',16);
    fact = c.factor_scores(:,2);
    if ~isempty(c.outcome)
        fact=[fact c.outcome]; fact=sortrows(fact,2);
        xlab='Component 2';
        
        lo = find(fact(:,2)<median(fact(:,2)));
        hi = find(fact(:,2)>median(fact(:,2)));
        
        plot(fact(lo,1),lo,'ko','LineWidth',2);
        plot(fact(hi,1),hi,'k^','LineWidth',2);
        set(gca,'YTick',[1 max(lo)+.5 size(fact,1)],'YTickLabel',{'Lowest' 'Median' 'Highest'})
        title(xlab)
        ylabel('Behavior');
        xlabel('Score')
    else
        xlab='Scores on Component 1';
        plot(fact(:,1),'ko','LineWidth',2);
        title(xlab); ylabel('Index'); xlabel('Scores: Comp. 2')
    end
    
end % if dosideplots

figure(f1)
axis off
drawnow

save_figure('cluster_nmdsfig',500);

return


% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Re-plot an nmds figure with outcome and covs added (but not part of
% clustering solution); Use existing completed c structure

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
function c_complete = nmdsfig_plot_withcovs(c, cl, doranks)

extravars = [c.outcome c.covs_nointerest];
nextra = size(extravars, 2);

% get data again, and do not remove covariates
c_complete = c;
c_complete = rmfield(c_complete, 'covs_nointerest');
c_complete = nmdsfig_tools('get_data_matrix',c_complete,cl,'timeseries',1,[],doranks);

disp('Creating c_complete.dat with [Outcome Covs Vars]');
c_complete.dat = [extravars c_complete.dat];

c_complete = nmdsfig_tools('get_correlations',c_complete);
[c_complete.GroupSpace,c_complete.obs,c_complete.implied_dissim,c_complete.stress] = shepardplot(c_complete.D,[]);

c_complete.ClusterSolution = c.ClusterSolution;

cl_id = max(c_complete.ClusterSolution.classes) + 1;  % class ID for extra points

c_complete.ClusterSolution.classes = [cl_id .* ones(nextra, 1); c_complete.ClusterSolution.classes];

extranames = {};
if ~isempty(c.outcome), extranames{1} = 'Outcome'; end
if ~isempty(c.covs_nointerest)
    for i = 1:size(c.covs_nointerest, 2)
        extranames{end + 1} = ['Cov' num2str(i)];
    end
end

c_complete.names = [extranames c.names];

nmdsfig_tools('nmdsfig_plot',c_complete, 0, 0, 'nofill');

saveas(gcf,'NMDS_complete_flat_plot_orig','png');

% direct effects, but should exclude OUTCOME?
[c_complete.direct_mtx, c_complete.mediated_mtx, c_complete.mediators] = matrix_direct_effects(c_complete.STATS.sigmat2, c_complete.dat);



nmdsfig_tools('removelines');
nmdsfig_tools('drawlines',c_complete.GroupSpace, c_complete.direct_mtx);
saveas(gcf,'NMDS_complete_flat_plot_direct','png');




return





% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Get correlation matrix, depending on ancova status

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


function c = get_correlations(c, varargin)

disp('Saving correlation matrix in c.r')

if ~isfield(c, 'ancova_codes'), c.ancova_codes = []; end

corrtype = 'r';
if length(varargin) > 1, corrtype = varargin{2}; end
if ~strcmp(corrtype, 'r') && ~isempty(c.ancova_codes)
    warning('nmdsfig_tools:incompat_inputs', 'Corr type cannot be entered with ANCOVA covariates.');
end

disp(['Correlation/association metric: ' corrtype])

if isempty(c.ancova_codes)
    doancova = 0;
    % no between covariates, just use correlations
    if strcmp(corrtype, 'r')
        [r,p]=corrcoef(c.dat); r(r>.9999) = 1;
    elseif strcmp(lower(corrtype), 'spearman')
        [r,p] = corr(c.dat, 'type', 'Spearman');
    else
        [r,t,p] = correlation(corrtype, c.dat);
    end
    
    c.STATS.corrtype = corrtype;
    
    [pthr,sig] = fdr_correct_pvals(p,r);
    
    c.STATS.p_thr = pthr;
    
else
    %doancova = input('Use covariate scores as covariate in ANCOVA? (inter-region correls then done w/i group): 1/0: ');
    doancova = 1;
end

if doancova
    
    disp('Using ANCOVA to model between-subjects effects')
    disp('Inter-region correlations within groups; interaction tests group x brain covariance interaction')
    % main effect of grp in this case is hi vs lo behavior on region Y (brain act.) after
    % controlling for region X -- a bit hard to interpret
    
    if length(unique(c.ancova_codes) > 2),
        disp('Binarizing ANCOVA code regressor into two groups.');
        c.ancova_codes = mediansplit(c.ancova_codes);
    end
    
    [b,t,p] = ancova(c.ancova_codes,[],c.dat); % first cell is grp (diff), 2nd is r, 3rd is r diff x group
    pall = p;   % save all slopes
    
    r = b{2};   % overall slope
    r = r + eye(size(r));   % add 1 to diagonal
    p = pall{2};
    [c.STATS.pthr(1),sig] = fdr_correct_pvals(p,r);
    
    ri = b{3};  % behavior x covariance interaction
    %ri = ri + eye(size(ri));   % add 1 to diagonal
    pi = pall{3};
    [c.STATS.pthr(2),sigi] = fdr_correct_pvals(pi,ri);
    
    c.STATS.siginteract = sigi;
    c.STATS.b = b; c.STATS.t = t; c.STATS.p = pall;
    c.STATS.btp_descrip = 'first cell is grp (diff), 2nd is r, 3rd is r diff x group.  betas, t-vals, p-vals';
    
else
    c.STATS.p = p;
end

c.doancova = doancova;
c.r = r;
c.STATS.sigmat = sig;
c.STATS.sigmat2 = (p < .05) .* sign(r);
c.STATS.p_descrip = 'sig. p-values are FDR corrected.';
c.STATS.p_descrip2 = 'sigmat2 is p < .05 uncorrected.';

fprintf('Positive connections (q < .05 FDR): %3.0f\n', sum(c.r(find(c.STATS.sigmat(:))) > 0))
fprintf('Negative connections (q < .05 FDR): %3.0f\n', sum(c.r(find(c.STATS.sigmat(:))) < 0))

disp('Getting dissimilarity matrix from correlation matrix, saving in c.D')
% scale correlations (or beta weights) so that they are a similarity matrix between 0 and 1.
c.D = (1 - c.r) ./ 2;
%c.D=distosim(c.D);
% Make sure it's a dissim. matrix
c.D=simtodis(c.D);

return




function [pthr,sig] = fdr_correct_pvals(p,r)

psq = p; psq(find(eye(size(p,1)))) = 0;
psq = squareform(psq);
pthr = FDR(p,.05);
if isempty(pthr), pthr = 0; end

sig = sign(r) .* (p < pthr);

return



% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Find cluster solution

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


function c = cluster_solution(c, objectsByDimsData, n_cluster_range, niter, names)
% -----------------------------------------------------------------------
% testcluster:
%
% Permute objects in space, get null hypothesis cluster quality, and test
% observed cluster solution against this.
% -----------------------------------------------------------------------

% pval is p-values for each number of clusters tested
% classes = class ("network") assignments for best clustering solution

c.ndims = size(objectsByDimsData, 2);
c.GroupSpace = objectsByDimsData;

[c.ClusterSolution.pvals, ...
    c.ClusterSolution.classes, ...
    c.ClusterSolution.classnames, ...
    c.ClusterSolution.X, ...
    c.ClusterSolution.include, ...
    c.ClusterSolution.names, ...
    c.ClusterSolution.stats] = ...
    testclustnew(objectsByDimsData, n_cluster_range, c.ndims, niter, names, 'keep', 'average');

c.ClusterSolution.niter = niter;

return




% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Apply clusuter solution and get class average data

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


function c = apply_clusters(c, varargin)

dointeractive = 1;
if length(varargin) > 0, dointeractive = varargin{1}; end
%
% requires c. ...
% ClusterSolution.classes
%
% uses c. ...
% outcome
%
% creates
% c.ClusterSolution.classnames
% c.class_avg_data

% add names of regions if missing
if ~isfield(c, 'names') || isempty(c.names) || length(c.names) < size(c.dat, 2)
    disp('Missing or invalid region names; adding.');
    for ii = 1:size(c.dat, 2)
        c.names{ii} = ['V' num2str(ii)];
    end
end

% name 'networks'
% -----------------------------------------------------------------------

if ~isfield(c,'ClusterSolution'), error('You need a ClusterSolution field first! ... try testclustnew'), end

if ~isfield(c,'dat'), error('Enter data for each variable in .dat field first'), end

if ~isfield(c,'colors') || isempty(c.colors)
    %c.colors = {'yo' 'bo' 'go' 'ro' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
    c.colors = {[1 1 0] [0 0 1] [0 1 0] [1 0 0] [1 .5 0] [0 1 1] [1 0 1] [.5 1 0] [0 .5 1] [0 1 .5] [1 0 .5]};
end

if ~isfield(c,'GroupSpace')
    disp(c)
    fieldnm = input('Type name of GroupSpace field with MCA/MDS/PCA coordinates: ', 's');
    myGroupSpace = c.(fieldnm);
else
    myGroupSpace = c.GroupSpace;
end

if ~isfield(c,'doancova')
    c.doancova = 0;
end

classnames = [];
if isfield(c,'APPLY_CLUSTER') && isfield(c.APPLY_CLUSTER,'names')
    classnames = c.APPLY_CLUSTER.names;
end

if (isempty(classnames) || length(classnames) < max(c.ClusterSolution.classes)) && (isempty(dointeractive) || dointeractive)
    donames = input('Name Clusters (Classes)? (1/0) ');
    
    if donames
        for i = 1:max(c.ClusterSolution.classes)
            disp(['Network ' num2str(i)])
            tmp = find(c.ClusterSolution.classes == i);
            for j = 1:length(tmp)
                fprintf(1,'%s . ',c.names{tmp(j)});
                if j > 4, fprintf(1,'\n'),end
            end
            fprintf(1,'\n')
            classnames{i} = input('Name this network: ','s');
        end
        
    else
        for i = 1:max(c.ClusterSolution.classes)
            classnames{i} = ['Net' num2str(i)];
        end
    end
end

% if we have a between-subjects covariate, enter it, otherwise empty
if ~isfield(c,'outcome'), c.outcome = []; end

% Apply this clustering solution to the clusters
% Get individual 'network' scores and use them as predictors of behavior
% in multiple regression
% -----------------------------------------------------------------------

% now run for the average across task conditions
% average correlations within and between networks, and test
if isfield(c, 'r')
    mycorrmtx = c.r; disp('Using .r field for correlation matrix.');
elseif isfield(c, 'corr')
    mycorrmtx = c.corr; disp('Using .corr field for correlation matrix.');
else
    error('Need .r or .corr field in input structure.')
end

c.APPLY_CLUSTER = apply_cluster_solution(c.ClusterSolution.classes,...
    mycorrmtx,...
    'names',classnames,'bcov',c.outcome, 'n', size(c.dat, 1), 'dointeractive', dointeractive);

c.APPLY_CLUSTER.names = classnames;

% matrix plot is not meaningful here; close
close;
%save_figure('cluster_nmdsfig_rmatrix',400);

% standardize variables?  No, let's let low variance vars contribute less
% to avgs

c.class_avg_dat = get_class_avgs(c,c.dat);
[c.APPLY_CLUSTER.avgr,c.APPLY_CLUSTER.avgp] = corrcoef(c.class_avg_dat);
c.APPLY_CLUSTER.avgsig = c.APPLY_CLUSTER.avgp < .055 .* sign(c.APPLY_CLUSTER.avgr);

% get class (network) average stimulus locations (group space) and plot
% -----------------------------------------------------------------------

Gs = myGroupSpace(:,1:c.ndims);
cla = c.ClusterSolution.classes;
for i = 1:max(cla)
    wh = find(cla == i);
    if length(wh)>1;
        classGs(i,:) = mean(Gs(wh,:));
    else
        classGs(i,:) = Gs(wh,:);
    end
end
c.APPLY_CLUSTER.Gs = classGs;
c.APPLY_CLUSTER.classes = 1:max(cla);

% significance of r among class averages (done above)
% [r,p]=corrcoef(c.class_avg_dat); r(r>.9999) = 1;
% sig = (p < .05) .* sign(r);
%[pthr,sig] = fdr_correct_pvals(p,r);

% Figures
if strcmp(c.doancova,'Yes'), create_wide_figure; subplot(1,2,1); else create_figure; end
set(gca,'FontSize',14)
nmdsfig(classGs(:,1:2),'names',c.APPLY_CLUSTER.names,'classes',c.APPLY_CLUSTER.classes, ...
    'sig', c.APPLY_CLUSTER.avgsig, ... %c.APPLY_CLUSTER.group_stats.sigu,'sizescale',[4 16]);
    'sizescale',[4 16],'colors',c.colors);

title(['Correlations among class averages (p <= .05 two-tailed)'])

if strcmp(c.doancova,'Yes'), subplot(1,2,2);
    set(gca,'FontSize',14)
    mdsfig(classGs(:,1:2),classnames,cla,c.STATS.siginteract);
    title(['Group by cov interaction: ANCOVA'])
end

save_figure('cluster_nmdsfig_class_lines_uncor',550);

return


% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Get class average data

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

function class_avg_dat = get_class_avgs(c,dat)

class_avg_dat = [];

[nobs,nvars] = size(dat);
nclasses = length(c.ClusterSolution.classes);
if nvars ~= nclasses, error('Number of classes does not match num. columns in dat'); end

for i = 1:max(c.ClusterSolution.classes)
    wh = find(c.ClusterSolution.classes == i);
    if length(wh) > 1
        class_avg_dat = [class_avg_dat nanmean(dat(:,wh)')'];
    else
        class_avg_dat = [class_avg_dat dat(:,wh)];
    end
end

return

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Stepwise regressions on behavior

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

function c = predict_behavior(c,meth)
% needs:
% c.
% outcome
% meth is 'classes' or 'regions'

if isempty(c.outcome)
    disp('Between-subjects c.outcome is empty.  No correlations with behavior run.');
    return
end

switch meth
    case 'classes'
        
        % Add logistic regression for categorical 1/0 DVs!
        % see also stepwise_tor.m which does this stuff.
        disp(' ')
        disp('------------------------------------------')
        disp('Using networks to predict behavior')
        disp('------------------------------------------')
        disp(' ')
        str = correlation_to_text([c.outcome c.class_avg_dat],[],[{'Behav'} c.APPLY_CLUSTER.class_names]);
        
        %[c.ClusterSolution.STEPWISE.b,c.ClusterSolution.STEPWISE.se,c.ClusterSolution.STEPWISE.pval, ...
        %c.ClusterSolution.STEPWISE.inmodel,stats] = ...
        %stepwisefit(c.class_avg_dat,c.outcome);
        disp(' ')
        c.class_STEPWISE = stepwise_tor(c.class_avg_dat,c.outcome,c.APPLY_CLUSTER.class_names);
        STEP = c.class_STEPWISE;
        names = c.APPLY_CLUSTER.class_names;
        
    case 'regions'
        
        %N = stats.dfe+stats.df0;
        %r2 = (stats.SStotal-stats.SSresid) ./ stats.SStotal;
        %adjr2 = 1 - (1-r2)*((N - 1) ./ (N - stats.df0 - 1));
        %fprintf(1,'Omnibus F(%3.0f,%3.0f) = %3.2f, RMSE = %3.2f, p = %3.4f, Adj R^2 = %3.2f\n', ...
        %stats.df0,stats.dfe,stats.fstat, ...
        %stats.rmse,stats.pval,adjr2);
        %stats.adjr2 = adjr2;
        %c.ClusterSolution.STEPWISE.stats = stats;
        
        disp(' ')
        disp('------------------------------------------')
        disp('Using individual regions to predict behavior'),
        disp('------------------------------------------')
        disp(' ')
        str = correlation_to_text([c.outcome c.dat],[],[{'Behav'} c.names],1);
        disp(' ')
        c.indiv_STEPWISE = stepwise_tor(c.dat,c.outcome,c.names);
        STEP = c.indiv_STEPWISE;
        names = c.names;
        
    otherwise error('Unknown method.  Valid = classes or regions')
end

% ------------------------------------------
% make plots
% ------------------------------------------
doplot = 1;

if ~isfield(c,'colors')
    c.colors = {'yo' 'bo' 'go' 'ro' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
end
while length(c.colors) < size(c.dat, 2), c.colors = [c.colors c.colors]; end

wh = logical(STEP.inmodel);

switch meth
    case 'classes'
        X = c.class_avg_dat(:, wh);
    case 'regions'
        X = c.dat(:, wh);
end


names = names(:, wh);

colors_to_use = c.colors(wh);

nsig = size(X,2);

if doplot && nsig > 0
    
    disp('Plotting partial correlation scatterplots for significant variables.');
    
    % Check for existing figure, and create if necessary
    create_figure('class_scatterplots', 1, nsig);
    
    scn_export_papersetup(300);
    
    for i = 1:nsig
        subplot(1, nsig, i);
        
        prplot(c.outcome, X, i, 0, colors_to_use(i));
        title(names{i})
        xlabel(' ');
        ylabel(' ');
        
    end
    
end


return



% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Draw lines on nmdsfig graph

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

function line_handles = drawlines(pc,sigmat,varargin)
% function handles_out = drawlines(pc,sigmat,[color],[style],[bendval])
%
% pc is stim coords
% sigmat is matrix of which lines to draw
% sigcol is multiplier to scale line width (0 is fixed width)
%   sigcol is a scalar between zero and one, lower is 'more salient'
% legmat is not used now (for legend stuff)
%
% Example:
% [hh3,hh4] = drawlines(pc,wh,sigcol(i),legmat,[.5 .5 .5; 0 1 1],[{':'} {':'}]);

% sigcol is a scalar between zero and one, lower is 'more salient'
sigcol = 0;
hold on

% linew
lw = 2;

color(1,:) = [0 0 0];   % first line color, 'positive'
color(2,:) = [0 1 1];   % second line color, 'negative'
style{1} = '-';         % first line style
style{2} = '-';         % second line style

if length(varargin) > 0 && ~isempty(varargin{1}), color = varargin{1}; end
if length(varargin) > 1 && ~isempty(varargin{2}), style = varargin{2}; sigcol = 0; end

bendval = .1;
if length(varargin) > 2 && ~isempty(varargin{3}), bendval = varargin{3};  end

hhp=[];hhn=[];
for i = 1 : size(pc,1)
    for j = i+1 : size(pc,1)
        
        if sigmat(i,j) > 0 || sigmat(i,j) < 0
            out = nmdsfig_tools('connect2d', [pc(i,1) pc(i,2)], [pc(j,1) pc(j,2)], 'b', lw - (lw*sigcol)-1, bendval);
            
            if sigmat(i,j) > 0
                hhp(end+1) = out.h;
                %hhp(end+1) = line([pc(i,1) pc(j,1)],[pc(i,2) pc(j,2)]);
                set(hhp,'Color',color(1,:) + [1 1 1] * abs(sigcol),'LineStyle',style{1},'LineWidth',lw - (lw*sigcol)-1)
            elseif sigmat(i,j) < 0
                hhn(end+1) = out.h;
                %hhn(end+1) = line([pc(i,1) pc(j,1)],[pc(i,2) pc(j,2)]);
                set(hhn,'Color',color(2,:) + abs([sigcol sigcol 0]),'LineStyle',style{2},'LineWidth',lw - (lw*sigcol)-1)
            end
            
        end
    end
end
%legend ([hhp hhn],legmat);

% store in gui data for later deletion, etc.
figh = gcf;
linehandles = [hhp hhn];
data = guidata(figh);
if ~isfield(data,'linehandles')
    data.linehandles = linehandles;
else
    data.linehandles = [data.linehandles linehandles];
end
guidata(figh,data);

line_handles.hhp = hhp;
line_handles.hhn = hhn;

return


% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Remove lines from nmdsfig plot

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
function doneok = removelines

doneok = 0;
figh = findobj('Tag', 'nmdsfig');
data = guidata(figh);
if isfield(data,'linehandles') && ~isempty(data.linehandles)
    delete(data.linehandles(ishandle(data.linehandles)))
    data.linehandles = [];
    guidata(figh,data);
    doneok = 1;
else
    disp('Cannot find line handles to remove in guidata.  Doing nothing.')
end

return



% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Put stars on figure where there are sig. relations in stepwise regression
% also do stars for robust correlations in glm (sep. function)
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
function sig = stepwise_stars(c,meth)

switch meth
    case 'regions'
        myfield = 'indiv_STEPWISE';
        myf2 = 'GroupSpace';
        
    case 'classes'
        myfield = 'class_STEPWISE';
        myf2 = 'APPLY_CLUSTER.Gs';
        
    otherwise
        error('Unknown method.  regions or classes.')
end

sig = c.(myfield).inmodel' .* c.(myfield).t;
plot_stars(c,sig,myf2);

return


function sig = robust_correl_stars(c,wh_reg)

[b,t,p,sig,F,fp,fsig,stat] = robust_reg_matrix([c.X ones(size(c.X,1),1)],c.dat,1);
sig = sig(wh_reg,:)' .* t(wh_reg,:)';
plot_stars(c,sig,'GroupSpace');

return


function plot_stars(c,sig,myfield)

wh = find(sig);
nsig = length(wh);

if nsig == 0
    disp('No significant results to plot.')
    return
end

% coordinates
pc = c.(myfield);
pc = pc(wh,1:2);
signval = (sig(wh) > 0) + 1;    % 1 for negative, 2 for positive
colors = {'bp' 'rp'};

for i = 1:nsig
    plot(pc(i,1),pc(i,2),colors{signval(i)},'MarkerSize',18,'LineWidth',3);
end

return


% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Get unique colors for each object, but similar colors within a class

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

function colors = class_colors(classes,basecolors)


nclasses = length(unique(classes));
nobj = length(classes);

if nargin < 2 || isempty(basecolors), basecolors = {[1 0 0] [0 1 0] [0 0 1] [1 1 0] [0 1 1] [1 0 1]}; end

while nclasses > length(basecolors), basecolors = [basecolors basecolors]; end

colors = cell(1,nobj);
mmax = .4;       % max for off-color value
%omin = .8;       % min for on-color value
rands = .25;

for i = 1:max(classes)
    
    wh = find(classes == i);
    
    n = length(wh);
    
    x = sort(linspace(0,mmax,n),'descend')';
    
    c = repmat(x,1,3);
    
    % random variation
    c = c + rands .* (rand(size(c)) .* sign((rand(size(c)) > .5) - .5));
    c(c<0) = 0;
    c(c>1) = 1;
    
    % add basecolor back in
    whone = find(basecolors{i} > 0);
    c(:,whone) = 1;
    
    for j = 1:n
        colors(wh(j)) = {c(j,:)};
    end
    
end

%     for i = 1:nobj
% fprintf(1,'%3.0f\t%3.2f %3.2f %3.2f\n',MCA.ClusterSolution.classes(i),colors{i})
% end
return



% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Utility functions

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% This is now its own separate function
% % function f1 = create_figure(tagname)
% %     % checks for old figure, and creates new one
% %
% %     if nargin < 1 || isempty(tagname)
% %         tagname = 'nmdsfig';
% %     end
% %
% %     old = findobj('Tag', tagname);
% %
% %     if ~isempty(old)
% %         disp(['Clearing old image: ' tagname]);
% %         clf(old);
% %     end
% %
% %     scnsize = get(0,'ScreenSize');
% %
% %     xdim = min(scnsize(3)./2, 700);
% %     ydim = min(scnsize(4)./2, 700);
% %
% %     f1 = figure('position',round([50 50 xdim ydim]),'color','white');
% %     set(f1, 'Tag', tagname);
% %
% %     set(gca, 'FontSize', 16);
% %
% %
% %     return


function f1 = create_wide_figure
scnsize = get(0,'ScreenSize');

myposition = [50 50 scnsize(3)-100 scnsize(4)/2];
myposition(3) = min(myposition(3), 1200);
myposition(4) = min(myposition(4), 500);

f1 = figure('position',myposition,'color','white');
return


function save_figure(myname,npts)
if nargin < 2, npts = 400; end
disp(['Saving figure: ' myname '.png'])
f1 = gcf;
scn_export_papersetup(npts);
saveas(f1,myname,'png');
return




% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Connect two points in 3-D space with a straight or curved line

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

function out = connect3d(varargin)

% ----------------------------------------------------------------------
% Parse inputs
% ----------------------------------------------------------------------

p = inputParser;

% Validation functions - customized for each type of input
% ----------------------------------------------------------------------

valfcn_scalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar'});

valfcn_number = @(x) validateattributes(x, {'numeric'}, {'nonempty'}); % scalar or vector

% Validation: Region object, atlas, structure, or [x1 x2 x3] triplet
valfcn_custom = @(x) isstruct(x) || isa(x, 'region') || isa(x, 'atlas') || (~isempty(x) && all(size(x) - [1 3] == 0) && all(isnumeric(x)));

% Validation: [x1 x2 x3] triplet
valfcn_xyz = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'size', [1 3]});

% Required inputs
% ----------------------------------------------------------------------
p.addRequired('x', valfcn_custom);
p.addRequired('y', valfcn_custom);

% Optional inputs
% ----------------------------------------------------------------------
% Pattern: keyword, value, validation function handle

p.addParameter('color', [.9 .2 0], valfcn_xyz);
p.addParameter('bendpercent', .15, valfcn_number); % can be scalar or vector
p.addParameter('thickness', .1, valfcn_scalar);
p.addParameter('nstreamlines', 30, valfcn_scalar);
p.addParameter('streamlineshift', 10, valfcn_scalar);
p.addParameter('nsamples', [], valfcn_scalar);

% Parse inputs and distribute out to variable names in workspace
% ----------------------------------------------------------------------
% e.g., p.parse([30 1 0], [-40 0 10], 'bendpercent', .1);
p.parse(varargin{:});

IN = p.Results;
fn = fieldnames(IN);

for i = 1:length(fn)
    str = sprintf('%s = IN.(''%s'');', fn{i}, fn{i});
    eval(str)
end

% Adjust input vars as needed
% ----------------------------------------------------------------------

if isa(x, 'atlas'), x = atlas2region(x); end
if isa(y, 'atlas'), y = atlas2region(y); end

% x, y inputs can be [x,y,z] coordinate triplets or clusters.
if isstruct(x) || isa(x, 'region'), x = x.mm_center; end
if isstruct(y) || isa(y, 'region'), y = y.mm_center; end

% make x, y, z bend percents
if length(bendpercent) == 1
    bendpercent = repmat(bendpercent,1,3);
end

% make 3 coords, so we can bend if we want to
xcoords = [x(1); x(1)+(y(1)-x(1))./2 + bendpercent(1)*x(1); y(1)];
ycoords = [x(2); x(2)+(y(2)-x(2))./2 + bendpercent(2)*x(2); y(2)];
zcoords = [x(3); x(3)+(y(3)-x(3))./2 + bendpercent(3)*x(3); y(3)];

n = length(xcoords);
 
if isempty(nsamples)
    nsamples = 10 * n;
end

% Plot streamlines
% ----------------------------------------------------------------------

for i = 1:nstreamlines % number of streamlines
    
    randval = unifrnd(-.15, .15, 4, 1);
    
    mycolor = color + randval(1);
    mycolor(mycolor > 1) = 1;
    mycolor(mycolor < 0) = 0;
    
    if any(bendpercent)
        % bow out: curved line
        
        myx = xcoords;
        myy = ycoords;
        myz = zcoords;
     
        % random bend
        myx(2) = myx(2) + streamlineshift*randval(2);
        myy(2) = myy(2) + streamlineshift*randval(2);
        myz(2) = myz(2) + streamlineshift*randval(2);
        
        t = 1:n;
        ts = 1:((n-1)/(nsamples-1)):n;          % spline grid
        
        myx = spline(t, myx, ts);
        myy = spline(t ,myy, ts);
        myz = spline(t, myz, ts);
    end
    
    % random shift
    myx = myx + streamlineshift * randval(2);
    myy = myy + streamlineshift * randval(3);
    myz = myz + streamlineshift * randval(4);
    
    h(i) = plot3(myx, myy, myz, 'Color', mycolor,'LineWidth', thickness);
    
end

% h = plot3(xcoords,ycoords,zcoords,'Color',color,'LineWidth',thickness);

out.h = h;
out.xcoords = xcoords;
out.ycoords = ycoords;
out.zcoords = zcoords;

return

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Connect two points in 2-D space with a straight or curved line

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
function out = connect2d(x,y,color,thickness,bendpercent,varargin)

% x, y inputs can be [x,y,z] coordinate triplets or clusters.
if isstruct(x), x = x.mm_center(1:2); end
if isstruct(y), y = y.mm_center(1:2); end

% make x, y bend percents
if length(bendpercent) == 1
    bendpercent = repmat(bendpercent,1,2);
end

% make 2 coords, so we can bend if we want to
xcoords = [x(1); x(1)+(y(1)-x(1))./2 + bendpercent(1)*x(1); y(1)];
ycoords = [x(2); x(2)+(y(2)-x(2))./2 + bendpercent(2)*x(2); y(2)];

if any(bendpercent)
    % bow out: curved line
    n = length(xcoords);
    
    nsamples = [];
    if length(varargin) > 0
        nsamples = varargin{1};
    end
    
    if isempty(nsamples)
        nsamples = 10 * n;
    end
    
    t = 1:n;
    ts = 1:((n-1)/(nsamples-1)):n;          % spline grid
    
    xcoords = spline(t,xcoords,ts);
    ycoords = spline(t,ycoords,ts);
    
end

h = plot(xcoords,ycoords,'Color',color,'LineWidth',thickness);

out.h = h;
out.xcoords = xcoords;
out.ycoords = ycoords;

return


