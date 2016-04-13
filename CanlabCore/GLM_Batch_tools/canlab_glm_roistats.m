function [ROISTATS] = canlab_glm_roistats(DIRS,ROI,varargin)
% Use canlab_glm_maskstats instead

fprintf('Use canlab_glm_maskstats instead\n');

return

% Returns a ROISTATS structure containing data from robfitdir's input subject level
% SPM analyses.
%
% :Usage:
% ::
%
%     ROISTATS = canlab_glm_roistats(robfitdir, roi)
%
% :Output structure:
%
%   **ROISTATS
%
%   **COV:**
%        subject x covariate matrix (EXPT.cov from robfit design)
%
%   **COVNAME:**
%        names of covariates in COV
%
%   **ROI:**
%        array of structs (1 per roi)
%
%   **MASK:**
%        the filename of the roi mask used
%
%   **SUB:**
%        struct containing data from subject level images
%
%   **CON:**
%        struct contains data from contrast images
%
%   **NAME:**
%        cell array (1 cell per contrast) of contrast names
%
%   **FILES:**
%        cell array (1 cell per contrast) of character arrays of 
%        contrast image filenames
%
%   **DATA:**
%        cell array (1 cell per contrast) of voxel x subject matrices of
%        contrast values (vectorized image data from dat field of fmri_data objects)
%
%   **MEANS:**
%        subject x contrast matrix of contrast means across voxels within
%        roi mask
%
%   **COVxCON:**
%        arrays of covariate matrix X contrast means correlation
%                   results (RHO and P, see help corr())
%
%   **GRP:**
%        struct containing data from group level images
%
%   **BETA:**
%        struct array (1 struct per group level regressor) of data
%        from beta images
%
%   **NAME:**
%        name of group level regressor
%
%   **FILES:**
%        cell array (1 cell per robust directory) of beta image files
%
%   **DATA:**
%        cell array (1 cell per robust directory) of vectorized image
%        data 
%
%   **MEANS:**
%        array (1 value per robust directory) of beta means 
%        across voxels within roi mask
%
% The following plots are saved in each roi's directory in the plots directory:
%   - contrast means by contrast (means are lines across subjects on x axis)
%   - contrast means by subject (means are dots, lined up along the x axis by contrast)
%   - group level betas (bar plot with group level regressors grouped by
%                        subject level contrasts)
%
% If there's more than one regressor in the group level model, for each regressor:
%   - scatter plot of subject level contrast means against group level regressor
%
% :Arguments:
%
%   **robfitdir:**
%        a directory in which robfit was run
%       contains robfit directories (e.g., robust0001)
%       preferably contains EXPT.mat or EXPTm.mat
%
%   **roi:**
%        a filename or cell array of filenames of masks to apply to data
%
% :Options:
%
%   **'c', connums:**
%        (vector of contrast numbers)
%       only include data from contrasts specified by numbers
%
%   **'c', conname(s):**
%        (string or cell array of strings)
%       only include data from contrasts specified by name
%
%   **'data':**
%        include DATA cells in output structure (they can take up space)
%
%   **'noplots':**
%        no plots will be made
%
%   **'od', dir:**
%        will save plots in dir (DEFAULT: robfitdir/roistats_scatterplots)
%
%   **'dotproduct':**
%        include dot product of weighted mask multiplication ("pattern expression")
%
% ..
%    Programmers' notes:
%    to write:
%    save data option
%    save as csv option
%    option to include subject-level beta data
%    option to break input mask into regions
%    within subjects error bars
%    grouped bar plots?
%    save volInfo (for mask? betas? cons? from fmri_data or spm?)
% ..


DO_PLOTS = true;
DO_SAVEDATA = false;
DO_DOTPRODUCT = false;

%% SETUP
% set defaults

applymaskopts = '';

% modify arguments as needed
if ~iscell(ROI), ROI = {ROI}; end

% error checking
if ~exist('DIRS','var') || isempty(DIRS)
    error('No robftidir specified.')
end
if ~exist('ROI','var') || isempty(ROI)
    error('No ROI(s) specified.')
end
for i=1:numel(ROI)
    if ~exist(ROI{i},'file')
        error('No such file: %s',ROI{i})
    end
end

% parse varargin
i=1;
while i<=numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'c'
                i=i+1;
                if ischar(varargin{i})
                    included_connames = {varargin{i}};
                elseif iscellstr(varargin{i})
                    included_connames = varargin{i};
                elseif isnumeric(varargin{i})
                    included_connums = varargin{i};
                else
                    error('argument to ''c'' must be number, string, or cell string')
                end
            case 'data'
                DO_SAVEDATA = true;
            case 'od'
                i=i+1;
                plotdir = varargin{i};
            case 'noplots'
                DO_PLOTS = false;
            case 'dotproduct'
                DO_DOTPRODUCT = true;
            case 'l1norm'
                applymaskopts = [applymaskopts ', ''l1norm'''];
            case 'l2norm'
                applymaskopts = [applymaskopts ', ''l2norm'''];
            otherwise
                error(['UNRECOGNIZED OPTION: ' varargin{i}])
        end
    else
        disp(varargin{i})
        error('Above option UNRECOGNIZED')
    end
    i=i+1;
end




%% PREP
if ischar(DIRS)
    DO_GRPLEV = true;
    BASEDOUT = DIRS;
    
    robdirs = filenames(fullfile(DIRS,'robust[0-9][0-9][0-9][0-9]'),'absolute');
    if isempty(robdirs)
        error('No such directories: %s/robust[0-9][0-9][0-9][0-9]',DIRS)
    end
    fprintf('... LOADING SETUP files from %d robust* directories in %s\n',numel(robdirs),DIRS)
    for r = 1:numel(robdirs)
        load(fullfile(robdirs{r},'SETUP'));
        SETUPS{r} = SETUP; %#ok
        allconnames{r} = SETUP.name;
        allconfiles{r} = SETUP.files;
        clear SETUP        
    end
    
    % load EXPT
    if exist(fullfile(DIRS,'EXPTm.mat'),'file')
        load(fullfile(DIRS,'EXPTm'));
        EXPT = EXPTm;
        clear EXPTm;
    elseif exist(fullfile(DIRS,'EXPT.mat'),'file')
        load(fullfile(DIRS,'EXPT'));
    end       
    
    % get covariates, if any
    if size(SETUPS{1}.X,2)>1
        fprintf('... GETTING between-subjects covariate(s)\n')
        cov = SETUPS{1}.X(:,2:end);
        
        covname = {};
        if isfield(EXPT,'covnames')
            covname = EXPT.covnames(2:end);
        end
    else
        cov = {};
        covname = '';
    end
    ROISTATS.cov = cov;
    ROISTATS.covname = covname;
elseif iscellstr(DIRS)
    DO_GRPLEV = false;
    BASEDOUT = pwd;
    
    for i = 1:numel(DIRS)    
        clear SPM spmmat
        spmmat = fullfile(DIRS{i},'SPM.mat');
        if isempty(spmmat)
            error('No such file: %s',spmmat)
        end
        load(spmmat);
        tmpconnames{i} = [SPM.xCon.name];
        allconfiles{i} = filenames(fullfile(DIRS{i}),'con_*.img');        
    end
    
    if size(unique(tmpconnames))~=1
        error('Not all subject levels have same contrast file names.')
    end    
    allconnames = cellstr(strvcat(SPM.xCon.name));
    clear SPM spmmat tmpconnames    
end


%% MAIN
for r=1:numel(ROI)
    clear mask
    
    mask = ROI{r};
    [ignore roiname] = fileparts(mask);
    fprintf('... GETTING data for roi: %s\n',mask)
    
    %% SUBJECT LEVEL STAT MAPS (cons)
    clear lconname lconfiles lcondata 
    
    %% load data      
    lcondotproduct = [];
    
    fprintf('Subject level contrasts:\n')
    n=0;
    if exist('included_connums','var')
        connums = included_connums;
    else
        connums = 1:numel(allconfiles);
    end
    for i = connums
        thisconname = deblank(allconnames{i});
        if exist('included_connames','var') && ~any(~cellfun('isempty',(regexp(thisconname,included_connames))))
            fprintf('EXCLUDED: %s\n',thisconname)
            continue
        else        
            n=n+1;
            included_robustdirs(n) = i; %#ok
        end
        % name
        lconname{n} = thisconname; %#ok
        fprintf('%s\n',lconname{n})
        % files
        lconfiles{n} = SETUPS{i}.files; %#ok
        % data
        clear dat
        evalc('dat = fmri_data(lconfiles{n},mask);');
        lcondata{n} = dat.dat; %#ok
        
        if DO_DOTPRODUCT % this needs to be done with the full fmri_data object, so this is the best place for it, though it belongs conceptually in the next section
            clear dotproduct
            evalc(['dotproduct = apply_mask(dat,mask,''pattern_expression''' applymaskopts ');']);
            lcondotproduct = [lcondotproduct dotproduct]; %#ok
        end
    end
    
    
    %% get subject level stats
    lconmeans = [];  
    fprintf('... GETTING stats from subject level contrasts\n')
    for i = 1:numel(lcondata)
        % means
        if size(lcondata{i},1) == 1
            lconmeans = [lconmeans lcondata{i}']; %#ok
        else
            lconmeans = [lconmeans nanmean(lcondata{i})']; %#ok
        end
    end    
            
    %% pack struct    
    ROISTATS.roi(r).mask = mask;
    ROISTATS.roi(r).sub.con.name = lconname;
    ROISTATS.roi(r).sub.con.files = lconfiles;
    ROISTATS.roi(r).sub.con.data = lcondata;
    ROISTATS.roi(r).sub.con.means = lconmeans;
    if DO_DOTPRODUCT, ROISTATS.roi(r).sub.con.dotproduct = lcondotproduct; end
    
    
    %% do t-testing on dotproduct values
    if DO_DOTPRODUCT
        [tt.hypothesis tt.p tt.ci tt.stats] = ttest(ROISTATS.roi(r).sub.con.dotproduct);
        ROISTATS.roi(r).sub.con.dotproduct_ttest = tt;
    end   
    
    %% make plots
    if DO_PLOTS        
        if ~exist('plotdir','var')  
            if DO_GRPLEV
                dout = fullfile(DIRS,'roistats_plots',roiname);
            else
                dout = fullfile(pwd,'roistats_plots',roiname);
            end
        else
            dout = fullfile(plotdir,roiname);
        end
        if ~exist(dout,'dir'), mkdir(dout); end
        
        fprintf('... MAKING plots\n')
        
        conticklabels = {};
        for i = 1:numel(ROISTATS.roi(r).sub.con.name)
            conticklabels{i} = sprintf('%s (%d)',ROISTATS.roi(r).sub.con.name{i},i); %#ok
        end
        
        subticklabels = {};
        for i = 1:size(ROISTATS.roi(r).sub.con.files{1},1)
            subticklabels{i} = regexprep(deblank(ROISTATS.roi(r).sub.con.files{1}(i,:)),'.*/([^/]*)/[^/]*','$1'); %#ok
            subticklabels{i} = sprintf('%s (%d)',subticklabels{i},i); %#ok
        end
        
        %%% means plots
        % contrast means across subjects
        figure; plot(ROISTATS.roi(r).sub.con.means);
        title(sprintf('subject level contrasts (mean within %s)',roiname),'interpreter','none');
        ylabel('beta contrast');        
        xaxis_labeling('subjects',subticklabels)
        fout = fullfile(dout,'contrast_means_by_subject.png');
        fprintf('saving: %s\n',fout);
        saveas(gcf,fout);
        close gcf
        
        % contrast means by contrast
        figure; plot(ROISTATS.roi(r).sub.con.means','o');
        title(sprintf('subject level contrasts (mean within %s)',roiname),'interpreter','none');
        ylabel('beta contrast');
        xaxis_labeling('contrasts',conticklabels)
        fout = fullfile(dout,'contrast_means_by_contrast.png');
        fprintf('saving: %s\n',fout);
        saveas(gcf,fout);
        close gcf
        
        if DO_DOTPRODUCT
            % dotproduct plot
            conticklabels = {};
            for i = 1:numel(ROISTATS.roi(r).sub.con.name)
                conticklabels{i} = sprintf('%s (%d)',ROISTATS.roi(r).sub.con.name{i},i); %#ok
            end
            
            figure; bar(nanmean(ROISTATS.roi(r).sub.con.dotproduct));
            title(sprintf('average dot product across subjects (weights: %s)',roiname),'interpreter','none');
            ylabel('dot product');
            xaxis_labeling('subject level contrasts',conticklabels);
            fout = fullfile(dout,'group_mean_dotproduct.png');
            fprintf('saving: %s\n',fout);
            saveas(gcf,fout);
            close gcf
        end
        
        
        %%% covariate plots
        if ~isempty(cov)
            [rho pval] = corr(cov,lconmeans);
            ROISTATS.roi(r).sub.con.covXcon.rho = rho;
            ROISTATS.roi(r).sub.con.covXcon.pval = pval;
            
            for i = 1:size(cov,2)
                for j = 1:size(lconmeans,2)
                    figure;
                    plot(cov(:,i),lconmeans(:,j),'o');
                    lsline;
                    xlabel(covname{i},'interpreter','none');
                    ylabel(sprintf('beta (meant within %s)',roiname),'interpreter','none');
                    title(sprintf('%s correlation with %s (r=%.4f, p<%.4f)',lconname{j},covname{i},rho(i,j),pval(i,j)),'interpreter','none');
                    
                    scon = regexprep(lconname{j},'[()]','');
                    scon = regexprep(scon,'[^A-Za-z0-9_+.-]','');
                    
                    gcov = regexprep(covname{i},'[()]','');
                    gcov = regexprep(gcov,'[^A-Za-z0-9_+.-]','');
                    
                    fout = fullfile(dout,[scon '__by__' gcov '.png']);
                    fprintf('saving: %s\n',fout)
                    saveas(gcf,fout);
                    close gcf
                end
            end
        end
    end
    
    
    %% GROUP LEVEL STAT MAPS (beta)
    
    %% load data
    clear hbetaname hbetafiles hbetadata hbetadataobj hbetameans hbetadotproduct
    
    hbetadotproduct = [];
    
    for b = 1:size(SETUPS{1}.X,2)
        % name
        if exist('EXPT','var') && isfield(EXPT,'covnames')
            hbetaname = EXPT.covnames{b};
        else
            hbetaname = 'hbeta1';
        end
                
        n=0;
        for i = included_robustdirs
            n=n+1;
            % file
            hbetafiles{n} = fullfile(robdirs{i},sprintf('rob_beta_%04d.img',b)); %#ok
            % data
            clear dat
            evalc('dat = fmri_data(hbetafiles{n},mask);');
            hbetadata{n} = dat.dat; %#ok
            
            if DO_DOTPRODUCT
                clear dotproduct
                evalc(['dotproduct = apply_mask(dat,mask,''pattern_expression''' applymaskopts ');']);
                hbetadotproduct = [hbetadotproduct dotproduct]; %#ok
            end
        end
        
        % means
        hbetameans = [];
        fprintf('... GETTING stats from group level beta: %s\n',hbetaname)
        for i = 1:numel(hbetadata)
            if size(hbetadata{i},1) == 1
                hbetameans = [hbetameans hbetadata{i}']; %#ok
            else
                hbetameans = [hbetameans nanmean(hbetadata{i})']; %#ok
            end
        end        
    
        % pack struct        
        ROISTATS.roi(r).grp.beta(b).name = hbetaname;
        ROISTATS.roi(r).grp.beta(b).files = hbetafiles;
        ROISTATS.roi(r).grp.beta(b).data = hbetadata;
        ROISTATS.roi(r).grp.beta(b).means = hbetameans;
        if DO_DOTPRODUCT, ROISTATS.roi(r).grp.beta(b).dotproduct = hbetadotproduct; end
    end
        
    if DO_PLOTS
        if ~exist('plotdir','var')
            dout = fullfile(BASEDOUT,'roistats_plots',roiname);
        else
            dout = fullfile(plotdir,roiname);
        end
        if ~exist(dout,'dir'), mkdir(dout); end
        
        fprintf('... MAKING plots\n')
        
        conticklabels = {};
        for i = 1:numel(ROISTATS.roi(r).sub.con.name)
            conticklabels{i} = sprintf('%s (%d)',ROISTATS.roi(r).sub.con.name{i},i); %#ok
        end
        
        X = [];
        for b = 1:size(SETUPS{1}.X,2)
            X = [X; ROISTATS.roi(r).grp.beta(b).means]; %#ok
        end
        X = X';                

        figure; bar(X,'group');
        title(sprintf('group level betas (mean within %s)',roiname),'interpreter','none');
        ylabel('beta');        
        legend(strvcat(ROISTATS.roi(r).grp.beta.name)); %#ok
        xaxis_labeling('subject level contrasts',conticklabels);
        fout = fullfile(dout,'group_level_betas.png');
        fprintf('saving: %s\n',fout);
        saveas(gcf,fout);
        close gcf
    end
    
    
    %% clean up
    close all        
end


%% save
if ~DO_SAVEDATA
    for r=1:numel(ROISTATS.roi)
        ROISTATS.roi(r).grp.beta = rmfield(ROISTATS.roi(r).grp.beta, 'data');
        ROISTATS.roi(r).sub.con = rmfield(ROISTATS.roi(r).sub.con, 'data');
    end
end

end


%% slightly edited code from mathworks.com (search "matlab angled tick labels")
function xaxis_labeling(xaxislabel,ticklabels)

% reduce size of axis to fit labels
pos = get(gca,'Position');
set(gca,'Position',[pos(1), .2, pos(3) .65])

% set tick locations
Xt = 1:numel(ticklabels);
Xl = [0 numel(ticklabels)+1];
set(gca,'XTick',Xt,'Xlim',Xl);


ax = axis; % Current axis limits
axis(axis); % Set the axis limit modes to manual
Yl = ax(3:4); % Y-axis limits

% Place the text labels
t = text(Xt,Yl(1)*ones(1,numel(Xt)),ticklabels,'Interpreter','none');
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
    'Rotation',45);

% Remove the default labels
set(gca,'XTickLabel','')

% get the Extent of each text object
for i = 1:length(ticklabels)
    ext(i,:) = get(t(i),'Extent'); %#ok
end

% Determine lower point for alignment of text
LowYPoint = min(ext(:,2));

% Place the axis label
XMidPoint = Xl(1) + abs(diff(Xl))/2;
text(XMidPoint,LowYPoint,xaxislabel,'VerticalAlignment','top','HorizontalAlignment','center');

% adjust width
OP = get(gca,'OuterPosition');
set(gcf,'OuterPosition',[OP(1) OP(2) numel(ticklabels)*50 OP(3)])

end
