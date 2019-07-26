function canlab_glm_subject_levels(dsgnarg, varargin)
% Performs lower level GLM analysis with SPM:
%   1. specifies model
%   2. estimates model
%   3. generates contrast images for model
%   4. creates directory with named links to spmT and con maps
%   5. publishes analyses with scn_spm_design_check
%
% :Usage:
% ::
%
%     canlab_glm_subject_levels(DSGN [options])
%
% DSGN struct - defines the model and analysis parameters
%
% canlab_glm_subject_levels('README') to see description
%
% :Options:
%
%   **'README':**
%        prints canlab_glm_README, an overview of canlab_glm_{subject,group}_levels
%
%   **'dsgninfo':**
%        prints description of DSGN structure
%
%   **'subjects', subject_list:**
%        ignore DSGN.subjects, use cell array subject_list
%
%   **'overwrite':**
%        turn on overwriting of existing analyses (DEFAULT: skip existing)
%
%   **'onlycons':**
%        only run contrast job (no model specification or estimation)
%        note: will overwrite existing contrasts
%        note: to not run contrasts, simply do not include a contrasts field in DSGN
%
%   **'addcons':**
%        only run contrasts that aren't already in SPM.mat
%        option to canlab_spm_contrast_job
%
%   **'nodelete':**
%        do not delete existing contrasts (consider using addcons, above)
%        option to canlab_spm_contrast_job
%
%   **'nolinks':**
%        will not make directory with named links to contrast images
%
%   **'noreview':**
%        will not run scn_spm_design_check
%
%   **'dream':**
%        if you're running on the dream cluster, this option will cause
%        all subjects to be run in parallel (submitted with matlab DCS and
%        the Sun Grid Engine)
%        Note: currently only works with MATLAB R2009a
%
%   **'email', address:**
%        send notification email to address when done running
%
% Model specification and estimation done by canlab_spm_fmri_model_job
%
% Contrasts are specified by canlab_spm_contrast_job_luka
% see that function for more info.
%
% ..
%     -------------------------------------------------------------------------
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
%
%     Programmers' notes:
% ..

STARTTIME = datestr(now,31);
STARTINGDIR = pwd;
%% SET UP

% set path
if exist('canlab_preproc.m','file') ~= 2
    try
        addpath(genpath('/data/projects/wagerlab/Repository'));
    catch %#ok
        error('Failed to add canlab repository (/data/projects/wagerlab/Repository) to path')
    end
end
addpath(genpath('/usr/local/spm/spm8/matlabbatch'));
addpath(genpath('/usr/local/spm/spm8/apriori'));
addpath(genpath('/usr/local/spm/spm8/canonical'));


% set defaults
diarydirname = 'canlab_glm_logs';
diaryfilename = ['subject_levels_' regexprep(regexprep(STARTTIME,' ','_'), ':','') '.log'];

OPTS.parallel_dream = '';
OPTS.overwrite = false;
OPTS.onlycons = false;
OPTS.nocatch = false;
OPTS.run_renaming = true;
OPTS.run_hist = false; % hidden option
OPTS.run_review = true;
OPTS.modeljob = '';
OPTS.conjob = '';


%% PARSE ARGUMENTS
% get DSGN structure
if ischar(dsgnarg)
    if strcmp(dsgnarg,'README')
        system(sprintf('cat %s',which('canlab_glm_README.txt')));
        return;
    elseif strcmp(dsgnarg,'dsgninfo')
        system(sprintf('cat %s',which('canlab_glm_dsgninfo.txt')));
        return;
    elseif any(regexp(dsgnarg,'\.mat$')) && exist(dsgnarg,'file')
        load(dsgnarg);
    else
        error('Unrecognized argument: %s', dsgnarg)
    end
elseif isstruct(dsgnarg)
    DSGN = dsgnarg;
else
    error('DSGN structure must be given as first argument, either as a variable or as a matfile.')
end

% go through varargin
i=1;
while i<=numel(varargin)
    if ischar(varargin{i})
        switch(varargin{i})
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
            case {'dsgninfo'}
                system(sprintf('cat %s',which('canlab_glm_dsgninfo.txt')));
                return;
            case {'subjects'}
                i=i+1;
                subjects = varargin{i};
            case {'overwrite'}
                OPTS.overwrite = true;
            case {'runhist'} % hidden option
                OPTS.run_hist = true;
            case {'onlycons'}
                OPTS.onlycons = true;
            case {'addcons'}
                OPTS.conjob = [OPTS.conjob ',''addcons'''];
            case {'nocatch'}
                OPTS.nocatch = true;
            case {'nodelete'}
                OPTS.conjob = [OPTS.conjob ',''nodelete'''];
            case {'nolinks'}
                OPTS.run_renaming = false;
            case {'noreview'}
                OPTS.run_review = false;
            case {'email'}
                i=i+1;
                address = varargin{i};
            otherwise
                error(['UNRECOGNIZED OPTION: ' varargin{i}])
        end
    else
        disp(varargin{i})
        error('Above option UNRECOGNIZED')
    end
    i=i+1;
end



%% DSGN-parsing
% catch bad fields
allowablefields = {...
    'metadata' ...
    'modeldir' ...
    'subjects' ...
    'funcnames' 'allowmissingfunc' ...
    'concatenation' ...
    'tr' 'hpf' 'fmri_t' 'fmri_t0' ...
    'conditions' 'pmods' 'convolution' 'multireg' 'singletrialsall' 'singletrials' 'ar1'...
    'allowmissingcondfiles' 'allowemptycond' 'notimemod' 'modelingfilesdir' ...
    'customrunintercepts' ...
    'contrasts' 'contrastnames' 'contrastweights' ...
    'regmatching' 'defaultsuffix' 'noscale' ...
    'timingcheck'
    };
actualfields = fieldnames(DSGN);
for i = 1:numel(actualfields)
    if isempty(strmatch(actualfields{i},allowablefields,'exact'))
        error('UNRECOGNIZED DSGN FIELD: DSGN.%s', actualfields{i})
    end
end


% parse inputs, set defaults for missing inputs, etc
if ~isfield(DSGN,'modeldir')
    error('No modeldir specified')
    %DSGN.modeldir = pwd;
end

if exist('subjects','var')
    DSGN.subjects = subjects;
else
    if ~isfield(DSGN,'subjects')
        error('No subjects specified')
    end
end
regexprep(DSGN.subjects,'/$','');

if ~isfield(DSGN,'contrastnames'), DSGN.contrastnames = {}; end
if ~isfield(DSGN,'contrastweights'), DSGN.contrastweights = {}; end

if isfield(DSGN,'defaultsuffix')
    OPTS.conjob = [OPTS.conjob ',''suffix'',''' DSGN.defaultsuffix ''''];
end

if isfield(DSGN,'noscale')
    if ~islogical(DSGN.noscale) && ~isnumeric(DSGN.noscale)
        error('DSGN.noscale must be true/false or 1/0')
    end
    if DSGN.noscale, OPTS.conjob = [OPTS.conjob ',''noscale''']; end
end

if isfield(DSGN,'regmatching')
    OPTS.conjob = [OPTS.conjob ',''' DSGN.regmatching ''''];
end


if ~OPTS.onlycons % check all the model spec/estimation stuff
    if ~isfield(DSGN,'funcnames'), error('No functional data specified'); end
    if ~isfield(DSGN,'allowmissingfunc'), DSGN.allowmissingfunc = false; end
    
    if isfield(DSGN,'timingcheck')
        if ~isfield(DSGN.timingcheck,'condition')
            error('Must define DSGN.timingcheck.condition');
        end
        
        if ~isfield(DSGN.timingcheck,'mask')
            error('Must define DSGN.timingcheck.mask');
        end
        if ~isfield(DSGN.timingcheck,'stat')
            DSGN.timingcheck.stat = 'mean';
        end
    end
        
    if ~isfield(DSGN,'conditions'), fprintf('WARNING: No conditions specified.\n'); end
    if ~isfield(DSGN,'allowmissingcondfiles'), DSGN.allowmissingcondfiles = false; end
    if ~isfield(DSGN,'allowemptycond'), DSGN.allowemptycond = false; end
    
    if ~isfield(DSGN,'convolution')
        DSGN.convolution.type = 'hrf';
        DSGN.convolution.time = 0;
        DSGN.convolution.dispersion = 0;        
    end        
    
    if ~isfield(DSGN,'singletrialsall'), DSGN.singletrialsall = false; end
    
    if ~isfield(DSGN,'ar1'), DSGN.ar1 = false; end
        
    if isfield(DSGN,'notimemod') && DSGN.notimemod
        OPTS.modeljob = [OPTS.modeljob ',''notimemod'''];
    end
    if ~isfield(DSGN,'modelingfilesdir'), DSGN.modelingfilesdir = 'spm_modeling'; end
    
    if ~isfield(DSGN,'tr'), error('no TR specified'); end
    if ~isfield(DSGN,'hpf'), error('no HPF specified'); end
    if isfield(DSGN,'fmri_t')
        OPTS.modeljob = [OPTS.modeljob ',''fmri_t'',' num2str(DSGN.fmri_t)];
    end
    if isfield(DSGN,'fmri_t0')
        OPTS.modeljob = [OPTS.modeljob ',''fmri_t0'',' num2str(DSGN.fmri_t0)];
    else
        error('DSGN.fmri_t0 (microtime onset) not specified')
    end
    
    if isfield(DSGN,'multireg')
        DSGN.multireg = [regexprep(DSGN.multireg, '\.mat$', '') '.mat'];
    else
        DSGN.multireg = '';
    end    
    
    if ~isfield(DSGN,'concatenation')
        DSGN.concatenation = {};
    end
    
    % check for uneven numbers of functional files
    for s = 1:numel(DSGN.subjects)
        for r = 1:numel(DSGN.funcnames)
            next(r) = numel(filenames(fullfile(DSGN.subjects{s},DSGN.funcnames{r}))); %#ok
        end
        if DSGN.allowmissingfunc
            if sum(next>1)~=0
                if ~isempty(DSGN.concatenation)
                    fprintf('ERROR: If allowmissingfunc is set and concatenation is desired,\n')
                    fprintf('  DSGN.funcnames strings must match one functional data file each\n')
                    return
                elseif numel(DSGN.conditions)>1
                    fprintf('ERROR: If allowmissing func is set:')
                    fprintf('  EITHER only one session''s worth of conditions is specified\n')
                    fprintf('  OR DSGN.funcnames strings match one functional data file each\n')
                    return
                end
            end
        elseif s>1 && any(next~=last)~=0
            fprintf('Number of functional images retrieved is not consistent across subjects.\n')
            fprintf('(fix this or use DSGN.allowmissingfunc)\n');
            return
        end
        last = next;
    end
end



%% PREP WORK
if ~exist(DSGN.modeldir,'dir')
    fprintf('Making subject level models directory: %s\n',DSGN.modeldir)
    mkdir(DSGN.modeldir)
else
    fprintf('Moving to existing subject level models directory: %s\n',DSGN.modeldir)
end
cd(DSGN.modeldir);

% initialize SPM's job manager
spm_jobman('initcfg');

% initialize statuses
modelstatus = zeros(1,numel(DSGN.subjects));
constatus = zeros(1,numel(DSGN.subjects));
linkstatus = zeros(1,numel(DSGN.subjects));
histstatus = []; % hidden option
reviewstatus = [];
tchkstatus = [];

% diary - to be on except when calling other functions that do their own logging
diarydir = fullfile(DSGN.modeldir,diarydirname);
if ~exist(diarydir,'dir'), mkdir(diarydir); end
diaryname = fullfile(diarydir,diaryfilename);
fulldiaryname = regexprep(diaryname,'\.log$','_full.log');
fprintf('Writing to logfile: %s\n',diaryname)
diary(diaryname), fprintf('STARTED: %s\n',STARTTIME), diary off

% make working directory
wd = regexprep(diaryname,'\.log$','');
mkdir(wd)



%% RUN LOWER LEVELS
if ~isempty(OPTS.parallel_dream)
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
    nargsout = 3;
    
    % submit jobs
    diary(diaryname), fprintf('%s\tSubmitting jobs to cluster\n',datestr(now,31)); date; diary off
    for s = 1:numel(DSGN.subjects)
        save(fullfile(wd,sprintf('env_%04d',s)),'DSGN','OPTS','STARTINGDIR');
        
        j = sched.createJob();
        set(j,'PathDependencies',cellstr(path));
        createTask(j, str2func('canlab_glm_subject_levels_run1subject'), nargsout, {wd s});
        alltasks = get(j, 'Tasks');
        set(alltasks, 'CaptureCommandWindowOutput', true);
        diary(diaryname), fprintf('Submitting analysis for %s\n',DSGN.subjects{s}); diary off
        submit(j);
    end
    
    % wait
    % look into using wait, waitForState
    diary(diaryname), fprintf('WAITING for jobs to stop running (they may run for a while)\n'); diary off        
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
    diary(diaryname), fprintf('\n%s\tJobs done running\n',datestr(now,31)); diary off
    
    % handle output arguments/streams
    diary(diaryname), fprintf('GATHERING job outputs\n'); diary off
    jobdirs = filenames(sprintf('%s/Job*[0-9]',wd),'absolute');
    for i = 1:numel(jobdirs)
        % load state and output
        jobin = load(sprintf('%s/Task1.in.mat',jobdirs{i}));
        jobout = load(sprintf('%s/Task1.out.mat',jobdirs{i}));
        
        % parse output arguments
        modelstatus(jobin.argsin{2}) = jobout.argsout{1};
        constatus(jobin.argsin{2}) = jobout.argsout{2};
        linkstatus(jobin.argsin{2}) = jobout.argsout{3};
        
        % output stream
        cw = jobout.commandwindowoutput;
        cw = regexprep(cw,'^> ','');
        cw = regexprep(cw,'\n> ','\n');
        cw = regexprep(cw,'\b[^\n]*\n',' LINES WITH BACKSPACES OMITTED\n');
        diary(fullfile(wd,sprintf('cmdwnd_%04d.txt',jobin.argsin{2}))), disp(cw), diary off
    end
    % output stream
    [ignore ignore] = system(sprintf('cat %s/cmdwnd_*txt > %s',wd,fulldiaryname)); %#ok
    
    % merge diaries
    [ignore ignore] = system(sprintf('rm %s/diary*',wd)); %#ok
    [ignore ignore] = system(sprintf('grep ''^> '' %s | sed ''s|^> ||'' >> %s',fulldiaryname,diaryname)); %#ok
else
    for s = 1:numel(DSGN.subjects)
        parsave(fullfile(wd,sprintf('env_%04d',s)),DSGN,OPTS,STARTINGDIR);
        [modelstatus(s) constatus(s) linkstatus(s)] = canlab_glm_subject_levels_run1subject(wd,s);
    end
    
    % merge diaries
    [ignore ignore] = system(sprintf('cat %s/diary* >> %s',wd,diaryname)); %#ok
end


%% POST ANALYSIS PROCESSES
diary(diaryname)
fprintf('\n\n\n')
fprintf('-------------------------------\n')
fprintf('--  POST ANALYSIS PROCESSES  --\n')
fprintf('-------------------------------\n')
diary off
if sum(modelstatus==-1) || sum(constatus==-1)
    diary(diaryname), fprintf('SKIPPED: some model jobs or contrast jobs failed.\n'), diary off
else
    % TIMING CHECK
    if isfield(DSGN,'timingcheck')
        diary(diaryname), fprintf('\n... GENERATING TIMING CHECK REPORTS\n'), diary off
        try           
            % publish_timing_check(DSGN);
            canlab_glm_subject_levels_timingcheck(DSGN);
            tchkstatus = 1;
        catch exc
            if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
            else diary(diaryname), fprintf('%s\n',getReport(exc,'extended')); diary off; end            
            tchkstatus = -1;
        end
    end
    
    if  sum(constatus==1) == 0
        diary(diaryname), fprintf('SKIPPED: no new contrasts have been run.\n'), diary off
    else
        % T MAP HISTOGRAMS
        if OPTS.run_hist % hidden option
            diary(diaryname), fprintf('\n... T STATISTIC HISTOGRAMS.\n'), diary off
            try
                cd(DSGN.modeldir)
                batch_t_histograms('o','t_histograms')
                close all
                histstatus = 1;
            catch exc
                if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
                else diary(diaryfile), fprintf('%s\n',getReport(exc,'extended')); diary off; end
                histstatus = -1;
            end
        end
        
        
        
        % CANLAB DESIGN REVIEW
        diary(diaryname), fprintf('\n... GENERATING DESIGN REVIEWS\n'), diary off
        if ~OPTS.run_review
            diary(diaryname), fprintf('SKIPPED: switched off\n'), diary off
        else
            try
                canlab_glm_publish('s',DSGN.modeldir);
                reviewstatus = 1;
            catch exc
                if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
                else diary(diaryfile), fprintf('%s\n',getReport(exc,'extended')); diary off; end
                reviewstatus = -1;
            end
        end
    end
end


%% failure report
diary(diaryname)
fprintf('\n\n\n')
fprintf('----------------------\n')
fprintf('--  FAILURE REPORT  --\n')
fprintf('----------------------\n')
nfail = failure_report(DSGN,modelstatus,constatus,linkstatus,histstatus,reviewstatus,tchkstatus);
diary off


%% email notification
if exist('address','var')
    try
        [ignore output] = system(sprintf('printf "canlab_glm_subject_levels has finished running.\nDirectory: %s\nLog file: %s\nFailure count: %d\n" | mail -v -s "canlab_glm_subject_levels done" %s',DSGN.modeldir,diaryname,nfail,address)); %#ok
    catch  %#ok
        diary(diaryname), fprintf('The notification email failed to send.\n'), diary off
    end
end


%% clean up
diary(diaryname), fprintf('\n\nFINISHED: %s\n',datestr(now,31)), diary off

cd(STARTINGDIR)

end




% -------------------------------------------------------------------------
%  SUBFUNCTIONS  ----------------------------------------------------------
% -------------------------------------------------------------------------

%%
function [nfail] = failure_report(DSGN,modelstatus,constatus,linkstatus,histstatus,reviewstatus,tchkstatus)

nfail = 0;
if sum(modelstatus == -1)
    fprintf('\nFAILED model spec/estim jobs:   %d of %d\n',sum(modelstatus==-1),sum(modelstatus~=0))
    fprintf('%s\n',DSGN.subjects{modelstatus == -1})
    nfail = nfail+sum(modelstatus == -1);
end
if sum(constatus == -1)
    fprintf('\nFAILED contrast jobs:   %d of %d\n',sum(constatus == -1),sum(constatus ~= 0))
    fprintf('%s\n',DSGN.subjects{constatus == -1})
    nfail = nfail + sum(constatus == -1);
end
if sum(linkstatus == -1)
    fprintf('\nFAILED named linking:   %d of %d\n',sum(linkstatus == -1),sum(linkstatus ~= 0))
    fprintf('%s\n',DSGN.subjects{linkstatus == -1})
    nfail = nfail + sum(linkstatus == -1);
end
if histstatus == -1, fprintf('\nFAILED t histogram making.\n'); nfail=nfail+1; end
if reviewstatus == -1, fprintf('\nFAILED design reviewing.\n'); nfail=nfail+1; end
if tchkstatus == -1, fprintf('\nFAILED timing check report.\n'); nfail=nfail+1; end

if ~nfail, fprintf('\nRAN WITH NO PROBLEMS (or at least so it seems).\n'); end

end


% function publish_timing_check(DSGN)
% 
% assignin('base','MASK',DSGN.timingcheck.mask);
% assignin('base','OP',DSGN.timingcheck.stat);
% 
% beforedir = pwd;
% 
% cd(DSGN.modeldir);
% 
% outputdir = fullfile(pwd,'timing_check');
% if exist(outputdir,'dir'), rmdir(outputdir,'s'); end
% mkdir(outputdir);
% 
% p = struct('useNewFigure', false, 'maxHeight', 1500, 'maxWidth', 1200, ...
%     'outputDir', outputdir, 'showCode', false);
% 
% fout = publish('canlab_glm_subject_levels_timingcheck.m',p);
% fprintf('Created subject level timing check:\n\t%s\n',fout);
% 
% cd(beforedir)
% 
% end

%% From Kathy Pearson:
% replace each first instance of SPM-like output backspace with newline;
% ignore additional backspaces found in sequence
%
function [wrapstr] = nobackspace(str) %#ok

wrapstr = str;
i = strfind(wrapstr, 8);
if ~isempty(i)
    k = 0;
    n = length(str);
    first8 = 1;
    for j = 1:n
        if str(j) == 8
            if first8
                k = k + 1;
                wrapstr(k) = 10;
                first8 = 0;
            end
        else
            k = k + 1;
            wrapstr(k) = str(j);
            first8 = 1;
        end
    end
    wrapstr = wrapstr(1:k);
end

end

function parsave(path,DSGN,OPTS,STARTINGDIR)
    save(path,'DSGN','OPTS','STARTINGDIR');
end