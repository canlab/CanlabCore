function [MASKSTATS] = canlab_glm_maskstats(DIRS,MASK,varargin)
% Returns a MASKSTATS structure containing data from robfitdir's input subject level
%    SPM analyses.
%
% :Usage:
% ::
%
%     MASKSTATS = canlab_glm_maskstats(robfitdir, mask, [options])
%
% output structure:
% MASKSTATS
%
%   **COV:**
%        subject x covariate matrix (EXPT.cov from robfit design)
%
%   **COVNAME:**
%        names of covariates in COV
%
%   **MASK:**
%        array of structs (1 per mask)
%
%   **MASKFILE:**
%        the filename of the mask used
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
%   **IMGFILES:**
%        cell array (1 cell per contrast) of character arrays of 
%        contrast image filenames
%
%   **(MEASURE):**
%        subject X contrast matrix of measures (see canlab_maskstats)
%
%   **COVxCON.(MEASURE):**
%        arrays of covariate matrix X contrast means correlation
%        results (RHO and P, see help corr())
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
%   **IMGFILES:**
%        cell array (1 cell per robust directory) of beta image files 
%
%   **(MEASURE):**
%        row vector (1 value per robust directory) of measures
%        (see canlab_maskstats)
%
% The following plots are saved in each mask's directory in the plots directory:
%   - contrast means by contrast (means are lines across subjects on x axis)
%   - contrast means by subject (means are dots, lined up along the x axis by contrast)
%   - group level betas (bar plot with group level regressors grouped by
%     subject level contrasts)
%
% If there's more than one regressor in the group level model, for each regressor:
%   - scatter plot of subject level contrast means against group level regressor
%
%
% :Arguments:
%
%   **robfitdir:**
%        a directory in which robfit was run
%        contains robfit directories (e.g., robust0001)
%        preferably contains EXPT.mat or EXPTm.mat
%
%   **mask:**
%        a filename or cell array of filenames of masks to apply to data
%
% :Options:
%
%   **MEASURE OPTIONS:**
%        (see canlab_maskstats) (DEFAULT: mean (within mask's non-zero voxels))
%
%   **'cons', connums:**
%        (vector of contrast numbers)
%        only include data from contrasts as specified by numbers
%
%   **'cons', conname(s):**
%        (string or cell array of strings)
%        only include data from contrasts as specified by name
%
%   **'plots':**
%        make plots
%
%   **'od', dir:**
%        will save plots in dir (DEFAULT: robfitdir/stats_scatterplots)
%
% ..
%    Programmers' notes:
%    to write:
%    save data option
%    save as csv option
%    option to include subject-level beta data
%    option to break input mask into regions?
%    within subjects error bars
%    grouped bar plots?
%    save volInfo (for mask? betas? cons? from fmri_data or spm?)
% ..


DO_PLOTS = false; % SETUP
% set defaults
OP = {};
ALLOPS = false;

applymaskopts = '';

% modify arguments as needed
if ~iscell(MASK), MASK = {MASK}; end

% error checking
if ~exist('DIRS','var') || isempty(DIRS)
    error('No robftidir specified.')
end
if ~exist('MASK','var') || isempty(MASK)
    error('No MASK(s) specified.')
end

% this gets checked/errored in canlab_maskstats
% for i=1:numel(MASK)
%     if ~exist(MASK{i},'file') && ~any(~cellfun('isempty',regexp({'nps','nps_thresh','nps_thresh_smooth'},['^' MASK{i} '$'])))
%         error('No such file: %s',MASK{i})
%     end
% end

% parse varargin
i=1;
while i<=numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'c' 'cons'}
                i=i+1;
                if ischar(varargin{i})
                    included_connames = varargin(i);
                elseif iscellstr(varargin{i})
                    included_connames = varargin{i};
                elseif isnumeric(varargin{i})
                    included_connums = varargin{i};
                else
                    error('argument to ''c'' must be number, string, or cell string')
                end
            case 'od'
                i=i+1;
                plotdir = varargin{i};
            case 'plots'
                DO_PLOTS = true;
            case {'binmean' 'mean' 'std' 'nonbinmean' 'nonbinstd' 'norm' ...
                  'dot_product' 'centered_dot_product' 'cosine_similarity' 'correlation'}
                OP{end+1} = varargin{i}; %#ok
            case 'all'
                ALLOPS = true;
            case 'l1norm'
                applymaskopts = [applymaskopts ', ''l1norm''']; %#ok
            case 'l2norm'
                applymaskopts = [applymaskopts ', ''l2norm''']; %#ok
            otherwise
                error(['UNRECOGNIZED OPTION: ' varargin{i}])
        end
    elseif iscellstr(varargin{i})
        for j=1:numel(varargin{i})
            switch varargin{i}{j}
                case {'binmean' 'mean' 'std' 'nonbinmean' 'nonbinstd' 'norm' ...
                        'dot_product' 'centered_dot_product' 'cosine_similarity' 'correlation'}
                    OP{end+1} = varargin{i}{j}; %#ok
                case 'all'
                    ALLOPS = true;
                otherwise
                    error('Unrecognized argument %s',varargin{i}{j})
            end
        end        
    else
        disp(varargin{i})
        error('Above option UNRECOGNIZED')
    end
    i=i+1;
end

if ALLOPS, OP = {'cosine_similarity' 'dot_product' 'correlation' 'mean' 'std' 'centered_dot_product'}; end
    
if isempty(OP), OP = {'mean'}; end




%% PREP
if ischar(DIRS)
    DO_GRPLEV = true;
    BASEDOUT = DIRS;
    
    robdirs = filenames(fullfile(DIRS,'robust[0-9][0-9][0-9][0-9]'),'absolute');
    if isempty(robdirs)
        error('No such directories: %s/robust[0-9][0-9][0-9][0-9]',DIRS)
    end
    fprintf('... LOADING SETUP files from %d robust* directories in %s\n',numel(robdirs),DIRS)
    SETUPS = cell(numel(robdirs),1); allconnames = SETUPS; allconfiles = SETUPS;
    for r = 1:numel(robdirs)
        load(fullfile(robdirs{r},'SETUP'));
        SETUPS{r} = SETUP;
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
    MASKSTATS.cov = cov;
    MASKSTATS.covname = covname;
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
        tmpconnames{i} = cellstr(char(SPM.xCon.name)); %#ok
        ntmpconnames(i) = size(tmpconnames{i},1); %#ok
        allconfiles{i} = filenames(fullfile(DIRS{i}),'con_*.img'); %#ok
    end
    
    if size(unique(ntmpconnames))~=1
        error('Not all subject levels have same contrast file names.')
    end    
    allconnames = tmpconnames{1};
    clear SPM spmmat tmpconnames    
end

if DO_PLOTS && ~exist('plotdir','var')
    plotdir = fullfile(BASEDOUT,'maskstats_plots');
end


%% MAIN
% for r=1:numel(MASK)
%     clear mask
%     
%     mask = MASK{r};
%     [ignore maskname] = fileparts(mask);
%     fprintf('... GETTING data for mask: %s\n',mask)
    
%% SUBJECT LEVEL STAT MAPS (cons)
fprintf('SUBJECT LEVELS\n')

clear lconname lconfiles lcondata

% get filenames
fprintf('Subject level contrasts:\n')
n=0;
lconfiles = {};
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
    lconfiles = [lconfiles cellstr(SETUPS{i}.files)]; %#ok
end

% get data
lmaskstats = canlab_maskstats(MASK,lconfiles,OP);

% pack output struct
for m=1:numel(MASK)
    MASKSTATS.mask(m).mask = lmaskstats(m).maskfile;
    MASKSTATS.mask(m).sub.con = lmaskstats(m).stats;
    MASKSTATS.mask(m).sub.con.imgfiles = lmaskstats(m).imgfiles;    
    MASKSTATS.mask(m).sub.con.name = lconname;
    if ~isempty(cov)
        for i=1:numel(OP)
            [rho pval] = corr(cov,MASKSTATS.mask(m).sub.con.(OP{i}));
            MASKSTATS.mask(m).sub.covXcon.(OP{i}).rho = rho;
            MASKSTATS.mask(m).sub.covXcon.(OP{i}).pval = pval;
        end
    end
end
    
    
% make plots
if DO_PLOTS
    fprintf('... MAKING plots\n')
    for m = 1:numel(MASK)
        [ignore maskname] = fileparts(MASK{m}); %#ok
        
        dout = fullfile(plotdir,maskname);
        if ~exist(dout,'dir'), mkdir(dout); end
        
        
        subticklabels = {};
        for i = 1:size(char(MASKSTATS.mask(m).sub.con.imgfiles{:,1}),1)
            subticklabels{i} = regexprep(deblank(MASKSTATS.mask(m).sub.con.imgfiles{i,1}),'.*/([^/]*)/[^/]*','$1'); %#ok
            subticklabels{i} = sprintf('%s (%d)',subticklabels{i},i); %#ok
        end
        
        for o=1:numel(OP)
            switch OP{o}
                case {'mean' 'dot_product' 'centered_dot_product' 'cosine_similarity' 'correlation'}
                    % line plot (one line for each contrast, across subjects)
                    
                    figure; plot(MASKSTATS.mask(m).sub.con.(OP{o}));
                    if strcmp(OP{o},'mean')
                        title(sprintf('subject level contrasts mean within %s',maskname),'interpreter','none');
                    else
                        title(sprintf('subject level contrasts %ss with %s',OP{o},maskname),'interpreter','none');
                    end
                    ylabel(sprintf('beta contrast %s',OP{o}),'interpreter','none');
                    xaxis_labeling('subjects',subticklabels)
                    fout = fullfile(dout,sprintf('contrasts_by_subject__%ss.png',OP{o}));
                    fprintf('saving: %s\n',fout);
                    saveas(gcf,fout);
                    close gcf
                    
                    % bar plot (mean across subjects per contrast)
                    evalc(['barplot_columns2(MASKSTATS.mask(m).sub.con.(OP{o}),maskname,'...
                        '''labels'',MASKSTATS.mask(m).sub.con.name,' ...
                        '''ylabel'',OP{o},''plabels'',''robust'',''ind'');']);
                    fout = fullfile(dout,sprintf('group_mean_%s.png',OP{o}));
                    fprintf('saving: %s\n',fout);
                    scn_export_papersetup(500);
                    saveas(gcf,fout);
                    close gcf
                    
                    %%% covariate plots
                    if ~isempty(cov)
                        subdout = fullfile(dout,'scatterplots');
                        if ~exist(subdout,'dir'), mkdir(subdout); end
                        for i = 1:size(cov,2)
                            for j = 1:size(MASKSTATS.mask(m).sub.con.(OP{o}),2)
                                figure;
                                plot(cov(:,i),MASKSTATS.mask(m).sub.con.(OP{o})(:,j),'o');
                                lsline;
                                xlabel(covname{i},'interpreter','none');
                                ylabel(sprintf('beta (meant within %s)',maskname),'interpreter','none');
                                title(sprintf('%s correlation with %s (r=%.4f, p<%.4f)',...
                                    MASKSTATS.mask(m).sub.con.name{j},covname{i},...
                                    MASKSTATS.mask(m).sub.covXcon.(OP{o}).rho(i,j),...
                                    MASKSTATS.mask(m).sub.covXcon.(OP{o}).pval(i,j)),...
                                    'interpreter','none');
                                
                                scon = regexprep(lconname{j},'[()]','');
                                scon = regexprep(scon,'[^A-Za-z0-9_+.-]','');
                                
                                gcov = regexprep(covname{i},'[()]','');
                                gcov = regexprep(gcov,'[^A-Za-z0-9_+.-]','');
                                
                                fout = fullfile(subdout,sprintf('%s_%ss__by__%s.png',scon,OP{o},gcov));
                                fprintf('saving: %s\n',fout)
                                saveas(gcf,fout);
                                close gcf
                            end
                        end
                    end
            end
        end
    end
end
    
    
%% GROUP LEVEL STAT MAPS (beta)
if DO_GRPLEV
    fprintf('GROUP LEVELS\n')
    nhcons = size(SETUPS{1}.X,2);
    hbetafiles = {};
    hbetaname = cell(nhcons,1);
    for b = 1:nhcons
        % get name
        if exist('EXPT','var') && isfield(EXPT,'covnames')
            hbetaname{b} = EXPT.covnames{b};
        else
            hbetaname{b} = sprintf('hbeta%04d',b);
        end
        
        % get list of files
        temphbetafiles = {};
        for i = included_robustdirs
            temphbetafiles{end+1,1} = fullfile(robdirs{i},sprintf('rob_beta_%04d.img',b)); %#ok
        end
        hbetafiles{b} = temphbetafiles; %#ok
    end
    hbetafiles = horzcat(hbetafiles{:});
    
    % extract data
    hmaskstats = canlab_maskstats(MASK,hbetafiles,OP);
        
    % pack output struct
    for b = 1:nhcons
        for m=1:numel(MASK)            
            for o=1:numel(OP)
                MASKSTATS.mask(m).grp.beta(b).(OP{o}) = hmaskstats(m).stats.(OP{o})(:,b);
            end
            MASKSTATS.mask(m).grp.beta(b).imgfiles = hmaskstats(m).imgfiles;
            MASKSTATS.mask(m).grp.beta(b).name = hbetaname{b};
        end
    end

    
    
    if DO_PLOTS
        dout = fullfile(plotdir,maskname);
        if ~exist(dout,'dir'), mkdir(dout); end
    
        fprintf('... MAKING plots\n')
        for m = 1:numel(MASK)            
            conticklabels = {};
            for i = 1:numel(MASKSTATS.mask(m).sub.con.name)
                conticklabels{i} = MASKSTATS.mask(m).sub.con.name{i}; %#ok
            end

            for o = 1:numel(OP)
                X = [];
                for b = 1:nhcons
                    X = [X MASKSTATS.mask(m).grp.beta(b).(OP{o})]; %#ok
                end
                
                figure; bar(X,'group');
                title(sprintf('group level betas %s with %s',OP{o},maskname),'Interpreter','none');
                ylabel(sprintf('beta %s',OP{o}),'Interpreter','none');
                legend(strvcat(MASKSTATS.mask(m).grp.beta.name)); %#ok
                xaxis_labeling('subject level contrasts',conticklabels);
                fout = fullfile(dout,sprintf('group_level_betas_%s.png',OP{o}));
                fprintf('saving: %s\n',fout);
                saveas(gcf,fout);
                close gcf
            end
        end
    end
end

%% clean up
close all

if DO_PLOTS
    % save MASKSTATS    
    for m = 1:numel(MASK)
        [ignore maskname] = fileparts(MASK{m}); %#ok
        dout = fullfile(plotdir,maskname);
        maskstats = MASKSTATS.mask(m); %#ok
        save(fullfile(dout,'raw_data'),'maskstats');
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
