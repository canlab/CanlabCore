function canlab_glm_group_levels(varargin)
% Performs group level robust GLM analysis with robfit
%   1. sets up analysis
%   2. runs robfit
%   3. (optionally) make inverse p maps (for FSL viewing)
%   4. (optionally) estimate significant cluster sizes
%   5. publishes analysis with robfit_results_batch
%
% :Usage:
% ::
%
%     canlab_glm_group_levels([options])
%
% :Optional Inputs:
%
%   **'s', subjects:**
%        cell array of filenames of subject-level analysis directories IN
%        modeldir
%        (note: modeldir won't be prepended to absolute paths)
%
%   **'m', modeldir:**
%        filename of directory containing subject level analyses
%
%   **'o', grpmodeldir:**
%        output directory name
%
%   **'c', cov:**
%        a matrix describing group level model
%        (do not include intercept, it is automatically included as first regressor)
%
%        see help robfit
%
%        note: requires specifying an output directory name
%
%        note: ordering of inputs (rows) must match subjects ordering
%
%   **'n', covname:**
%        a cell array of names for the covariates in cov
%
%   **'f', covfile:**
%        a csv file will specify the group level model:
%        first column, header 'subject', contains names of subject
%        directories (to be found in modeldir)
%        subsequent columns have covariates, headers are names of covariates
%        name of covfile will be name of group analysis directory (placed in
%        grpmodeldir)
%
%   **'mask', maskimage:**
%        filename of mask image
%
%   **DSGN:**
%        will use the following fields of the DSGN structure:
%        modeldir    = DSGN.modeldir
%
%        subjects    = DSGN.subjects
%
%        maskimage   = DSGN.mask
%
% :Note:
%   A covfile will cause other specifications of subject, cov, and 
%   covnames to be ignored.
%
%   If parameters are defined more than once (e.g., modeldir or subjects), 
%   only the last entered option will count.
%
% :Defaults:
%
%   **subjects:**
%        all SPM.mat-containing directories in modeldir  
%
%   **modeldir:**
%        pwd
%
%   **grpmodeldir:**
%        modeldir/one_sample_t_test
%
%   **cov:**
%        {} (run 1 sample t-test, see help robfit)
%
%   **covname:**
%        'groupmean'
%
%   **mask:**
%        'brainmask.nii'
%
% :Options:
%
%   **'README':**
%        prints canlab_glm_README, an overview of canlab_glm_{subject,group}_levels
%
%   **'overwrite':**
%        overwrite existing output directories
%
%   **'noresults':**
%        don't run/publish robfit_results_batch
%
%   **'onlyresults':**
%        just run/publish robfit_results_batch, don't run robfit (assumes existing analyses)
%
%   **'whichcons', [which cons]
%       vector of contrasts to analyze (DEFAULT: aall subject level contrasts)
%       see [which cons] in help robfit
%
%   **'invp' [, target_space_image]:**
%       generate inverse p maps and resample to the voxel and image dimensions
%       of target_space_image
%       (viewable on dream with ~ruzicl/scripts/invpview)
%
%   **'nolinks':**
%        do not make directory of named links to robust directories (using contrast names)
%
%   **'dream':**
%        if you're running on the dream cluster, this option will cause 
%       all analyses (e.g., lower level contrasts) to be run in parallel 
%       (submitted with matlab DCS and the Sun Grid Engine)
%       Note: currently only works with MATLAB R2009a
%
%   **'email', address:**
%       send notification email to address when done running
%
% ..
%     ----------------------------------------------------------------------
%     Copyright (C) 2013  Luka Ruzic
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
%     ----------------------------------------------------------------------
%
%     Programmers' notes
%     NOT READY YET:
%     'grf'
%       estimate significant cluster sizes for each contrast using GRF theory
%       will run at a=.05, voxelwise thresholds of .05, .01, .005, .001
%       see help estimate_cluster_extent
% ..

STARTTIME = datestr(now,31); % SET UP

STARTINGDIR = pwd;

% path
addpath(genpath('/usr/local/spm/spm8/matlabbatch'));
addpath(genpath('/usr/local/spm/spm8/apriori'));
addpath(genpath('/usr/local/spm/spm8/canonical'));

% set defaults
OPTS.parallel_dream = '';
OPTS.overwrite = false;
OPTS.run_results = true;
OPTS.run_robfit = true;
OPTS.run_printconnames = true;
OPTS.run_invp = false;
OPTS.run_grf = false;
OPTS.run_links = true;
OPTS.nocatch = false;

pthresh = [.05 .01 .005 .001]; % for significant cluster extent estimation

EXPT.cov = [];
EXPT.mask = which('brainmask.nii');

% set up statuses
robfitstatus = 0;
invpstatus = 0;
grfstatus = 0;
linkdirstatus = 0;
reviewstatus = 0;

diarydirname = 'canlab_glm_logs';
diaryfilename = ['group_levels_' regexprep(regexprep(STARTTIME,' ','_'), ':','') '.log'];

%% parse arguments
i=1;
while i<=numel(varargin)
    if isstruct(varargin{i})
        DSGN = varargin{i};
        if isfield(DSGN,'modeldir')
            modeldir = DSGN.modeldir;
        end
        if isfield(DSGN,'subjects')
            for s = 1:numel(DSGN.subjects)
                [dummy subname{s}] = fileparts(DSGN.subjects{s}); %#ok
                sublevs{s} = subname{s}; %#ok
            end
        end
        if isfield(DSGN,'mask')
            EXPT.mask = DSGN.mask;
        end
    elseif ischar(varargin{i})
        switch varargin{i}
            case {'README'}
                system(sprintf('cat %s',which('canlab_glm_README.txt')));
                return;
            case {'dream'}
                version = ver('matlab');
                switch version.Release
                    case '(R2009a)'
                        OPTS.parallel_dream = '2009a';
                    case '(R2011a)'
                        OPTS.parallel_dream = '2011a';
                    otherwise
                        error('Current version of matlab (%s)',version.Release)
                end
            case {'f'}
                i=i+1;
                covfile = varargin{i};
            case {'c'}
                i=i+1;
                EXPT.cov = varargin{i};
            case {'n'}
                i=i+1;
                EXPT.covnames = ['intercept' varargin{i}];
            case {'s'}
                i=i+1;
                sublevs = varargin{i};
            case {'d' 'm'}
                i=i+1;
                modeldir = varargin{i};
            case {'mask'}
                i=i+1;
                EXPT.mask = varargin{i};
            case {'overwrite'}
                OPTS.overwrite = true;
            case {'noresults'}
                OPTS.run_results = false;
            case {'onlyresults'}
                OPTS.run_robfit = false;
            case {'whichcons'}
                i=i+1;
                includedcons = varargin{i};
            case {'invp'}
                OPTS.run_invp = true;
            case {'grf'}
                OPTS.run_grf = true;
            case {'nolinks'}
                OPTS.run_links = false;
            case {'nocatch'}
                OPTS.nocatch = true;
            case {'o'}
                i=i+1;
                grpmodeldir = varargin{i};
            case {'email'}
                i=i+1;
                address = varargin{i};
            otherwise
                error(['Unrecognized argument: ' varargin{i}])
        end
    else
        disp(varargin{i})
        error('Above argument is unrecognized.')
    end
    i=i+1;
end

% MODELDIR
if ~exist('modeldir','var')
    modeldir = pwd;
end
if ~exist(modeldir,'dir')
    error(['No such directory: ' modeldir])
else
    % make absolute (unsure if necessary)
    modeldir = filenames(modeldir,'absolute','char');
end


% SUBLEVS
if exist('sublevs','var')
    for s = 1:numel(sublevs)
        if ~strcmp(sublevs{s}(1),'/')
            sublevs{s} = fullfile(modeldir,sublevs{s}); %#ok
        end
    end
end


% GRPMODELDIR
if ~exist('grpmodeldir','var')
    if exist('modeldir','var')
        grpmodeldir = fullfile(modeldir,'one_sample_t_test');
    else
        grpmodeldir = fullfile(pwd,'one_sample_t_test');
    end
end
% make grpmodeldir
if ~exist(grpmodeldir,'dir')
    mkdir(grpmodeldir)
end
% make absolute 
grpmodeldir = filenames(grpmodeldir,'absolute','char');


% COV
if exist('covfile','var')
    if ~exist(covfile,'file')
        error(['No such file: ' covfile])
    else
        [d f e] = fileparts(covfile); %#ok

        % grpmodeldir
        grpmodeldir = fullfile(fileparts(grpmodeldir),f);        
        
        % sublevs
        tempsubsfile = fullfile(d,[f '_tempsubs.txt']);
        system(['cut -d, -f1 < ' covfile ' | tail -n +2 > ' tempsubsfile]);        
        sublevs = importdata(tempsubsfile); delete(tempsubsfile);
        if isnumeric(sublevs), sublevs = textscan(num2str(sublevs'), '%s'); sublevs=sublevs{1}; end            
        for i=1:numel(sublevs)
            sublevs{i} = fullfile(modeldir,sublevs{i});
        end                
        
        % cov, covnames
        tempcovsfile = fullfile(d,[f '_tempcovs.txt']);
        system(['cut -s -d, -f2- < ' covfile ' > ' tempcovsfile]);
        s = importdata(tempcovsfile); delete(tempcovsfile);
        if isfield(s,'data')
            EXPT.cov = s.data;
            EXPT.covnames = ['intercept' s.colheaders];
        else
            EXPT.covnames = {'groupmean'};
        end
    end
else
    if isempty(EXPT.cov)
        EXPT.covnames = {'groupmean'};
    else
        if ~exist('grpmodeldir','var')
            error('cov has been specified: must specify grpmodeldir.')
        end
        if isfield(EXPT,'covnames')
            if numel(EXPT.cov) ~= numel(EXPT.covnames)
                error('Number of covariate names (%d) does not match number of covariates (%d).',numel(EXPT.covnames),numel(EXPT.cov))
            end
        else
            EXPT.covnames{1} = 'intercept';
            for i=2:numel(EXPT.cov)
                EXPT.covnames{i} = sprintf('grpcov%d',i-1);
            end
        end
    end
end


% set up diary - to be on except when calling other functions that do their own logging
diarydir = fullfile(grpmodeldir,diarydirname);
if ~exist(diarydir,'dir'), mkdir(diarydir); end
diarydir = filenames(diarydir,'absolute','char');
diaryname = fullfile(diarydir,diaryfilename);
fulldiaryname = regexprep(diaryname,'\.log$','_full.log');
fprintf('Writing to logfile: %s\n',diaryname)
diary(diaryname), fprintf('> STARTED: %s\n',STARTTIME), diary off

% make working directory
wd = regexprep(diaryname,'\.log$','');
mkdir(wd)

%% SET UP AND RUN ROBFIT
diary(diaryname), announce_string('ROBUST REGRESSION'), diary off
if ~OPTS.run_robfit
    diary(diaryname), fprintf('> SKIPPED: turned off in options\n'); diary off
else
    % subject level analyses + 1 SPM.mat (assumed representative)
    if ~exist('sublevs','var')
        diary(diaryname), fprintf('> No subject levels specified: looking for all subject level analyses in:\n\t%s\n',modeldir); diary off
        spmfiles = filenames(fullfile(modeldir,'*/SPM.mat'),'absolute');
        if isempty(spmfiles)
            diary(diaryname), fprintf('> ERROR: No subject level analyses found.\n'), diary off
            cd(STARTINGDIR)
            return
        end
        for i=1:numel(spmfiles), [sublevs{i}] = fileparts(spmfiles{i}); end %#ok
        sublevs = unique(sublevs);
        spmfile = spmfiles{1};
    else
        spmfile = fullfile(sublevs{1},'SPM.mat');
    end
    
    diary(diaryname)
    fprintf('> Subject level analyses:\n')
    for i=1:numel(sublevs)
        fprintf('>  %-8d%s\n',i,sublevs{i})
    end
    diary off
    
    if exist(spmfile,'file')
        load(spmfile)
    else
        diary(diaryname), fprintf('> ERROR: No such file: %s\n',spmfile); diary off
        cd(STARTINGDIR)
        return
    end
    
    % connames, connums, P
    diary(diaryname), fprintf('> Loading contrast names from:\n> \t%s\n',spmfile); diary off
    tcons = strmatch('T', char(SPM.xCon.STAT));
    EXPT.SNPM.connames = strvcat(SPM.xCon(tcons).name); %#ok
    
    ncons = size(EXPT.SNPM.connames,1);
    EXPT.SNPM.connums = 1:ncons;
    EXPT.SNPM.P = {};
    for c = 1:ncons
        confiles = {};
        for s = 1:numel(sublevs)
            confilename = fullfile(sublevs{s}, sprintf('con_%04d.img', tcons(c)));
            confile = filenames(confilename, 'char', 'absolute');
            if exist(confile,'file')
                confiles{s} = confile; %#ok
            else
                diary(diaryname), fprintf('> ERROR: No such contrast file: %s\n',confilename); diary off
                cd(STARTINGDIR)
                return
            end
        end
        EXPT.SNPM.P{c} = strvcat(confiles); %#ok
    end
    
    
    %% RUN ROBFIT
    diary(diaryname), fprintf('> Group level analysis directory:\n> \t%s\n',grpmodeldir); diary off
    cd(grpmodeldir);
    
    save EXPT EXPT
    
    % which contrasts should be run
    if ~exist('includedcons','var')
        includedcons = [];
        for c = 1:numel(EXPT.SNPM.connums)
            if OPTS.overwrite || numel(filenames(fullfile(pwd,sprintf('robust%04d', c),'rob_p_*img'))) < size(EXPT.cov,2)+1
                includedcons = [includedcons, c]; %#ok
            end
        end
    end
    
    if isempty(includedcons)
        diary(diaryname), fprintf('> No robust regressions to run.\n'); diary off
    else
        diary(diaryname)
        fprintf('> Subject level contrasts:\n');
        for i=1:size(EXPT.SNPM.connames,1)
            fprintf('> ');
            if ~sum(find(includedcons == i)), fprintf('*EXCLUDED*  '); end
            fprintf(' %-4d%s\n',i,EXPT.SNPM.connames(i,:));
        end
        fprintf('> Group level predictors:\n');
        for i=1:numel(EXPT.covnames)
            fprintf('>  %-4d%s\n',i,EXPT.covnames{i});
        end
        diary off
        
        njobs = numel(includedcons) * (1+numel(EXPT.covnames));
        if ~isempty(OPTS.parallel_dream) && njobs>1  
            % define scheduler
            sched = findResource('scheduler', 'type', 'generic');
            set(sched,'ClusterSize', 1);
            set(sched,'DataLocation',wd);
            switch OPTS.parallel_dream
                case '2009a'
                    set(sched,'ClusterMatlabRoot','/usr/local/matlab/R2009a');
                    set(sched,'ParallelSubmitFcn',@sgeParallelSubmitFcn);
                    set(sched,'SubmitFcn',@sgeSubmitFcn);
                    set(sched,'DestroyJobFcn', @sgeDestroyJobFcn);
                    set(sched,'DestroyTaskFcn',@sgeDestroyTaskFcn);
                case '2011a'
                    set(sched,'ClusterMatlabRoot','/usr/local/matlab/R2011a');
                    set(sched,'ParallelSubmitFcn',@parallelSubmitFcn);
                    set(sched,'SubmitFcn',@distributedSubmitFcn);
                    set(sched,'DestroyJobFcn', @destroyJobFcn);
                    set(sched,'DestroyTaskFcn',@destroyTaskFcn);
            end                     
            set(sched,'ClusterOsType','unix');
            set(sched,'HasSharedFileSystem', true);
            nargsout = 1;
            
            % submit jobs
            diary(diaryname), fprintf('> %s\tSubmitting jobs to cluster\n',datestr(now,31)); date; diary off            
            for c = 1:numel(includedcons)
                %                 save(fullfile(wd,sprintf('env_%04d',c)),'EXPT','includedcons','OPTS','grpmodeldir','pthresh');
                save(fullfile(wd,sprintf('env_%04d',c)),'EXPT','includedcons','OPTS','grpmodeldir','STARTINGDIR');
                j = sched.createJob();
                set(j,'PathDependencies',cellstr(path));
                createTask(j, str2func('canlab_glm_group_levels_run1input'), nargsout, {wd c});
                alltasks = get(j, 'Tasks');
                set(alltasks, 'CaptureCommandWindowOutput', true);
                diary(diaryname), fprintf('> Submitting analysis for contrast number %d\n',includedcons(c)); diary off
                submit(j);
            end
            
            % wait
            diary(diaryname), fprintf('> WAITING for jobs to stop running (they may run for a while)\n'); diary off
            t1 = clock;
            notdone = true;
            fprintf('Elapsed minutes:         ')
            while notdone
                notdone = false;
                for i=1:8, fprintf('\b'); end; fprintf('%8.1f',etime(clock,t1)/60);
                pause(10)
                jobstatefiles = filenames(sprintf('%s/Job*.state.mat',wd),'absolute');
                for s = 1:numel(jobstatefiles)
                    jobstate = textread(jobstatefiles{s},'%s');
                    if strcmp(jobstate,'running') || strcmp(jobstate,'queued'), notdone = true; break; end
                end
            end
            diary(diaryname), fprintf('> \n> %s\tJobs done running\n',datestr(now,31)); diary off
            
            % handle output arguments/streams
            diary(diaryname), fprintf('> GATHERING job outputs\n'); diary off
            jobdirs = filenames(sprintf('%s/Job*[0-9]',wd),'absolute');
            for i = 1:numel(jobdirs)
                % load state and output
                jobin = load(sprintf('%s/Task1.in.mat',jobdirs{i}));
                jobout = load(sprintf('%s/Task1.out.mat',jobdirs{i}));
                
                % parse output arguments
                robfitstatuses(jobin.argsin{2}) = jobout.argsout{1}; %#ok
%                 grfstatuses(jobin.argsin{2}) = jobout.argsout{2}; %#ok
                                
                % output stream
                diary(fullfile(wd,sprintf('cmdwnd_%04d.txt',jobin.argsin{2}))), regexprep(jobout.commandwindowoutput,'\b[^\n]*\n',' LINES WITH BACKSPACES OMITTED\n'), diary off
            end
            % merge statuses
            if any(robfitstatuses==-1), robfitstatus = -1; else robfitstatus = 1; end
%             if any(grfstatuses==-1), grfstatus = -1; else grfstatus = 1; end
            
            % output stream
            [ignore ignore] = system(sprintf('cat %s/cmdwnd_*txt > %s',wd,fulldiaryname)); %#ok
            
            % merge diaries
            [ignore ignore] = system(sprintf('grep ''^> '' %s >> %s',fulldiaryname,diaryname)); %#ok
        else
            cmd = 'EXPT = robfit(EXPT, includedcons, 0, EXPT.mask);';
            diary(diaryname), disp(cmd); diary off
            try
                eval(cmd)
                robfitstatus = 1;
            catch exc
                if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc);
                else diary(diaryname), fprintf('> %s\n',getReport(exc,'extended')), diary off; end
                robfitstatus = -1;
            end
        end
    end
    
    
    %% post processes    
    % write contrast names to textfiles
    if OPTS.run_printconnames
        subconnamefile = fullfile(grpmodeldir,'subconnames.txt');
        grpconnamefile = fullfile(grpmodeldir,'grpconnames.txt');
        diary(diaryname), announce_string('WRITING CONTRAST NAMES FILES'), diary off
        diary(diaryname), fprintf('> %s\n> %s\n',subconnamefile,grpconnamefile), diary off
        
        fid = fopen(subconnamefile,'w');
        for i = 1:size(EXPT.SNPM.connames,1)
            fprintf(fid,'> %d %s\n',i,deblank(EXPT.SNPM.connames(i,:)));
        end
        fclose(fid);
        
        fid = fopen(grpconnamefile,'w');
        for i = 1:numel(EXPT.covnames)
            fprintf(fid,'> %d %s\n',i,EXPT.covnames{i});
        end
        fclose(fid);
    end         
    
    
    % estimate significant cluster sizes
    %     if ~OPTS.parallel_dream && OPTS.run_grf
    if OPTS.run_grf
        diary(diaryname), announce_string('SIGNIFICANT CLUSTER SIZE ESTIMATION using GRF'), diary off
        graymattermask = which('scalped_avg152T1_graymatter_smoothed.img'); %#ok used in eval
        try
            robdirs = filenames(fullfile(grpmodeldir,'robust[0-9][0-9][0-9][0-9]'),'absolute');
            for i = 1:numel(robdirs)
                load(fullfile(robdirs{i},'SETUP.mat'));
                
                cmd = 'sigclext = estimate_cluster_extent(.05, pthresh, SETUP.files, ''mask'', graymattermask);';
                diary(diaryname), fprintf('> %s\n',cmd), diary off
                eval(cmd);
                close all
                
                fout = fullfile(robdirs{i},'significant_cluster_extents_grf.txt');
                dlmwrite(fout,[pthresh' sigclext(:,1)],'precision','%g','delimiter',' '); %#ok
            end
            
            grfstatus = 1;
        catch exc
            if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
            else diary(diaryname), fprintf('> %s\n',getReport(exc,'extended')), diary off; end
            grfstatus = -1;
        end
    end
    
    
    % make invp maps
    if OPTS.run_invp
        cmd = sprintf('robust_results_make_invp_maps(''d'',''%s'');',grpmodeldir);
        diary(diaryname), announce_string('MAKING INVERSE P MAPS'), diary off
        try
            eval(cmd);
            invpstatus = 1;
        catch exc
            if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc);
            else diary(diaryname), fprintf('> %s\n',getReport(exc,'extended')), diary off; end
            invpstatus = -1;
        end
    end
    
    
    % make named links to robust directories
    if OPTS.run_links
        diary(diaryname), announce_string('MAKING NAMED LINKS TO ROBUST DIRECTORIES'), diary off
        try
            % get grpmodeldir without /../
            if ismac
                [dummy absgrpmodeldir] = system(sprintf('cd %s; pwd -P',grpmodeldir)); %#ok ~ doesn't work in older matlabs
            else
                [dummy absgrpmodeldir] = system(sprintf('readlink -f %s',grpmodeldir)); %#ok ~ doesn't work in older matlabs
            end
            absgrpmodeldir = deblank(absgrpmodeldir);
            % make output directory
            linkdir = fullfile(absgrpmodeldir,'named_robust_directories');
            mkdir(linkdir);
            % make links
            for i = 1:size(EXPT.SNPM.connames,1)
                % get absolute path of robust directory (with no /../)
                robdir = fullfile(absgrpmodeldir, sprintf('robust%04d',i));
                
                % clean up contrast name for filesystem
                linkname = deblank(EXPT.SNPM.connames(i,:));
                linkname = regexprep(linkname,'[()]','');
                linkname = regexprep(linkname,'[^0-9A-Za-z-_]','_');
                linkname = fullfile(linkdir,linkname);
                
                % make link               
                diary(diaryname)
                fprintf('> %20s: %s\n','robust directory',robdir)
                fprintf('> %20s: %s\n','named link',linkname)
                eval(sprintf('!ln -s %s %s',robdir,linkname))
                diary off
            end
            
            linkdirstatus = 1;
        catch exc
            if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc);
            else diary(diaryname), fprintf('> %s\n',getReport(exc,'extended')), diary off; end
            linkdirstatus = -1;            
        end
    end
            
    
end


%% PUBLISH RESULTS
diary(diaryname), announce_string('PUBLISH RESULTS'), diary off
if ~OPTS.run_results
    diary(diaryname), fprintf('> SKIPPED: turned off in options.\n'); diary off
elseif robfitstatus == -1
    diary(diaryname), fprintf('> SKIPPED: robfit failed.\n'); diary off
else
    try
        canlab_glm_publish('g', grpmodeldir);
        reviewstatus = 1;
    catch exc
        if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc);
        else diary(diaryname), fprintf('> %s\n',getReport(exc,'extended')), diary off; end
        reviewstatus = -1;
    end
end

%% failure report
diary(diaryname)
announce_string('FAILURE REPORT')
fail = 0;
if robfitstatus == -1
    fprintf('> \n> WARNING: robfit failed.\n')
    fail=fail+1;
end
if invpstatus == -1
    fprintf('> \n> WARNING: inverse p map making failed.\n')
    fail=fail+1;
end
if grfstatus == -1
    fprintf('> \n> WARNING: cluster size estimation failed.\n')
    fail=fail+1;
end
if linkdirstatus == -1
    fprintf('> \n> WARNING: named linking to robust directories failed.\n')
    fail=fail+1;
end
if reviewstatus == -1
    fprintf('> \n> WARNING: review publishing failed.\n')
    fail=fail+1;
end

if fail==0, fprintf('> \n> RAN WITH NO PROBLEMS (or at least so it seems).\n'); end
diary off


%% email notification
if exist('address','var')
    diary(diaryname) 
    try
        [ignore output] = system(sprintf('printf "canlab_glm_group_levels has finished running.\nDirectory: %s\n\nLog file: %s\nFailed parts: %d\n" | mail -v -s "canlab_glm_group_levels done" %s',grpmodeldir,diaryname,fail,address)); %#ok
    catch %#ok
        fprintf('The email notification job failed.\n');
    end
    diary off
end

%% clean up
diary(diaryname), fprintf('> \n> \n> FINISHED: %s\n',datestr(now,31)), diary off

cd(STARTINGDIR)

end



% -------------------------------------------------------------------------
%  SUBFUNCTIONS  ----------------------------------------------------------
% -------------------------------------------------------------------------


function announce_string(string)

s = sprintf('--  %s  --',string);
l = regexprep(s,'.','-');
fprintf('> \n> \n> \n> %s\n> %s\n> %s\n> \n',l,s,l);

end


%% From Kathy Pearson:
% replace each first instance of SPM-like output backspace with newline;
% ignore additional backspaces found in sequence
%
% function [wrapstr] = nobackspace(str)
% 
% wrapstr = str;
% i = strfind(wrapstr, 8);
% if ~isempty(i)
%     k = 0;
%     n = length(str);
%     first8 = 1;
%     for j = 1:n
%         if str(j) == 8
%             if first8
%                 k = k + 1;
%                 wrapstr(k) = 10;
%                 first8 = 0;
%             end
%         else
%             k = k + 1;
%             wrapstr(k) = str(j);
%             first8 = 1;
%         end
%     end
%     wrapstr = wrapstr(1:k);
% end
% 
% end
