function [modelstatus constatus linkstatus] = canlab_glm_subject_levels_run1subject(wd, s)
% child process of canlab_glm_subject_levels
% (see canlab_glm_README.txt for an overview)
%
% ..
%    Copyright (C) 2013  Luka Ruzic
% ..

load(fullfile(wd,sprintf('env_%04d',s)));
%% PREP

diaryfile = fullfile(wd,sprintf('diary_%04d.log',s));

batchname = 'spm_specify_and_estimate_model';

if ~isempty(OPTS.parallel_dream)  %#ok
    z = '> ';
else
    z = '';
end

% initialize statuses
modelstatus = 0;
constatus = 0;
linkstatus = 0;

%subject-specific setup
[ignore subnum] = fileparts(DSGN.subjects{s}); %#ok
submodeldir = fullfile(DSGN.modeldir, subnum);    % make it later when you know it's not getting skipped
batchfile = fullfile(submodeldir, batchname);

diary(diaryfile)
fprintf('%s\n%s\n%s\n',z,z,z);
fprintf('%s------------------------------\n',z);
fprintf('%s--  SUBJECT LEVEL ANALYSIS  --\n',z);
fprintf('%s------------------------------\n',z);
fprintf('%s\n%sOutput Directory:\n%s\t%s\n',z,z,z,submodeldir);
diary off


%% MODEL SPECIFICATION AND ESTIMATION JOBS
diary(diaryfile), fprintf('%s\n%s... MODEL SPECIFICATION AND ESTIMATION JOBS\n',z,z), diary off
if OPTS.onlycons
    diary(diaryfile), fprintf('%sSKIPPED: turned off in options.\n',z), diary off
else
    run_this_model = true;
    %{
    if exist(submodeldir,'dir')
        %if ~numel(filenames(fullfile(submodeldir,'beta_*.img')))
        try status = importdata(fullfile(submodeldir,'.ssglm_model_status')); catch, status = ''; end %#ok
        if strcmp(status,'started')
            diary(diaryfile), fprintf('%sDELETING existing analysis directory: unfinished.\n',z), diary off
            rmdir(submodeldir,'s')
        elseif numel(filenames(fullfile(submodeldir,'beta_*.img'))) == 0
            diary(diaryfile), fprintf('%sDELETING existing analysis directory: no betas.\n',z), diary off
            rmdir(submodeldir,'s')
        else
            if OPTS.overwrite
                diary(diaryfile), fprintf('%sOVERWRITING: analysis directory exists.\n',z), diary off
                rmdir(submodeldir,'s')
            else
                diary(diaryfile), fprintf('%sSKIPPED: analysis directory exists.\n',z), diary off
                run_this_model = false;
            end
        end
    end
    %}
    
    if run_this_model
        if ~exist(DSGN.subjects{s},'dir')
            diary(diaryfile), fprintf(sprintf('%sERROR: no such data directory: %s\n',z,DSGN.subjects{s})), diary off
            modelstatus = -1;
            return
        end
        mkdir(submodeldir);
        eval(sprintf('!echo started > %s',fullfile(submodeldir,'.ssglm_model_status')))
        save(fullfile(submodeldir,'DSGN'),'DSGN');
        
        %% GET RUNS
        diary(diaryfile), fprintf('%sADDING input functional data\n',z), diary off
        clear runs runs3d
        try
            diary(diaryfile)
            [runs runs3d] = find_runs(DSGN,s,z);
            diary off
        catch exc
            if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
            else diary(diaryfile), fprintf('> %s\n',getReport(exc,'extended')); diary off; end
            modelstatus = -1;
            return
        end
                
        diary(diaryfile)
        fprintf('Functional data:\n');
        for r = 1:numel(runs), fprintf('%s\t%3d\t%s\n',z,r,runs{r}); end
        diary off
        
        
        %% PARSE CONDITIONS, REGRESSORS        
        diary(diaryfile), fprintf('%sADDING conditions and regressors\n',z), diary off
        clear names onsets durations pmods multipleregressors
        try
            [names onsets durations pmods multipleregressors] = parse_conditions(DSGN,runs,z);
        catch exc
            if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
            else diary(diaryfile), fprintf('> %s\n',getReport(exc,'extended')); diary off; end
            modelstatus = -1;
            return
        end
        
        %% CONCATENATION (if desired)
        if ~isempty(DSGN.concatenation)
            diary(diaryfile), fprintf('%sCONCATENATING data according to DSGN.concatenation:\n',z), diary off
            try          
                diary(diaryfile)
                [runs3d names onsets durations pmods multipleregressors] = concatdata(DSGN,submodeldir,runs,runs3d,names,onsets,durations,pmods,multipleregressors,z); %#ok
                diary off
            catch exc
                if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
                else diary(diaryfile), fprintf('> %s\n',getReport(exc,'extended')); diary off; end
                modelstatus = -1;
                return
            end
        end
        
        %% CONVERT CONDITIONS TO SINGLE TRIALS ANALYSIS (if desired)
        if DSGN.singletrialsall || isfield(DSGN,'singletrials')
            diary(diaryfile), fprintf('%sCONVERTING conditions to single trials analysis (one condition per trial)',z), diary off
            try
                diary(diaryfile)
                [names onsets durations pmods] = convert_to_single_trials(DSGN,names,onsets,durations,pmods);
                diary off
            catch exc
                if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
                else diary(diaryfile), fprintf('> %s\n',getReport(exc,'extended')); diary off; end
                modelstatus = -1;
                return
            end
        end
        
        %% RUN SPECIFICATION AND ESTIMATION BATCH
        % convert into flat arrays + conditions_by_run for canlab_spm_fmri_model_job
        clear conditions_by_run
        try
            diary(diaryfile)            
            [runs runs3d names onsets durations pmods conditions_by_run OPTS] = prep_for_canlab_spm_fmri_model_job(DSGN,OPTS,runs,runs3d,names,onsets,durations,pmods,z); %#ok
            diary off
        catch exc
            if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
            else diary(diaryfile), fprintf('> %s\n',getReport(exc,'extended')); diary off; end
            modelstatus = -1;
            return
        end
        
        modeljobcmd = ['matlabbatch = canlab_spm_fmri_model_job(submodeldir, DSGN.tr, DSGN.hpf, runs3d, conditions_by_run, onsets, durations, names, multipleregressors, ''pmod'', pmods ' OPTS.modeljob ');'];
        diary(diaryfile), fprintf('%s\nRUNNING model\n\t%s\n',modeljobcmd,z), diary off
        try
            eval(modeljobcmd);
            save(batchfile, 'matlabbatch');
            spm_jobman('run', matlabbatch);
            eval(sprintf('!echo finished > %s',fullfile(submodeldir,'.ssglm_model_status')))
            modelstatus = 1;
        catch exc
            if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
            else diary(diaryfile), fprintf('> %s\n',getReport(exc,'extended')); diary off; end
            modelstatus = -1;
            return
        end    
    end
    close all
end



%% CONTRAST JOB
diary(diaryfile), fprintf('%s\n%s... CONTRAST JOB.\n',z,z), diary off
try status = importdata(fullfile(submodeldir,'.ssglm_contrast_status')); catch, status = ''; end %#ok
if ~isfield(DSGN,'contrasts') || ~numel(DSGN.contrasts)~=0
    diary(diaryfile), fprintf('%sSKIPPED: no contrasts specified.\n',z), diary off
else
    if modelstatus ~= 1 % otherwise the model ran and must try to run cons
        if modelstatus == -1
            diary(diaryfile), fprintf('%sSKIPPED: model failed.\n',z), diary off
            return
        elseif ~OPTS.onlycons && numel(filenames(fullfile(submodeldir, 'spmT_*.img'))) == numel(DSGN.contrasts) && ~strcmp(status,'started')
            diary(diaryfile), fprintf('%sSKIPPED: already previously run, (and onlycons is not set).\n',z), diary off
            return
        end
    end
    
    conjobcmd = ['canlab_spm_contrast_job_luka(submodeldir, DSGN.contrasts, ''names'', DSGN.contrastnames, ''weights'', DSGN.contrastweights' OPTS.conjob ');'];
    diary(diaryfile), disp(conjobcmd), diary off
    try
        eval(sprintf('!echo started > %s',fullfile(submodeldir,'.ssglm_contrast_status')))
        eval(conjobcmd);
        eval(sprintf('!echo finished > %s',fullfile(submodeldir,'.ssglm_contrast_status')))
        constatus = 1;
    catch exc
        if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
        else diary(diaryfile), fprintf('> %s\n',getReport(exc,'extended')); diary off; end
        constatus = -1;
        return
    end
end



%% MAKE NAMED LINKS TO T MAPS
if OPTS.run_renaming && constatus==1
    diary(diaryfile), fprintf('%s\n%s... MAKING named links to spmT maps.\n',z,z), diary off
    try
        link_to_stat_maps(submodeldir)
        linkstatus = 1;
    catch exc
        if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
        else diary(diaryfile), fprintf('> %s\n',getReport(exc,'extended')); diary off; end
        linkstatus = -1;
    end
end


end


% -------------------------------------------------------------------------
%  SUBFUNCTIONS  ----------------------------------------------------------
% -------------------------------------------------------------------------


%%
function [runs runs3d] = find_runs(DSGN,session,z)

r=1;
for f = 1:numel(DSGN.funcnames)
    % find runs using path from funcnames (wildcards optional)
    runstoadd = filenames(fullfile(DSGN.subjects{session}, DSGN.funcnames{f}));
    if isempty(runstoadd)
        if DSGN.allowmissingfunc
            fprintf('%sWARNING: no runs found with: %s\n',z,DSGN.funcnames{f})
            runs{r} = ''; %#ok
            runs3d{r} = ''; %#ok
            r=r+1;
            continue
        else
            error('%sno runs found with: %s\n\t(to allow, set DSGN.allowmissingfunc)',z,DSGN.funcnames{f})
        end
    end
    
    % expand
    for j = 1:numel(runstoadd)
        runs{r} = runstoadd{j}; %#ok
        runs3d{r} = cellstr(expand_4d_filenames(runs{r})); %#ok
        r=r+1;
    end
end

end


%%
function [names onsets durations pmods multipleregressors] = parse_conditions(DSGN,runs,z)

newpmod = struct('name', [], 'param', [], 'poly', []);
emptypmod = newpmod([]);

for session = 1:numel(runs)
    % skip empty runs
    if isempty(runs{session}), continue; end
    
    % if only first session described, use for all sessions
    if size(DSGN.conditions,2) == 1
        sess = 1;
    else
        sess = session;
    end
    
    % specify modeling files directory
    [funcdir] = fileparts(runs{session});
    mfdir = fullfile(funcdir, DSGN.modelingfilesdir);
    
    c = 1; % keeps track of condition number
    
    for i = 1:numel(DSGN.conditions{sess})
        clear matfiles        
        
        if ~isempty(regexp(DSGN.conditions{sess}{i},'[[]{}*?]','once'))
            wc = fullfile(mfdir,DSGN.conditions{sess}{i});
            matfiles = filenames(wc);
            if isempty(matfiles)
                if DSGN.allowmissingcondfiles
                    fprintf('%sWARNING: No conditions files found with wildcard: %s\n',z,wc)                    
                    continue
                else
                    error('%sNo conditions files found with wildcard: %s\n    (See DSGN.allowmissingcondfiles)',z,wc)
                end
            end
        else
            matfiles{1} = fullfile(mfdir,DSGN.conditions{sess}{i});
            matfiles{1} = [regexprep(matfiles{1},'.mat$','') '.mat'];
            if ~exist(matfiles{1},'file')
                if DSGN.allowmissingcondfiles
                    fprintf('%sWARNING: No such conditions file: %s\n',z,matfiles{1})
                    continue
                else
                    error('%sNo such conditions file: %s\n    (See DSGN.allowmissingcondfiles)',z,matfiles{1})
                end
            end
        end
            
        % get pmods
        if isfield(DSGN,'pmods')
            try
                currpmods = emptypmod;
                for j = 1:numel(DSGN.pmods{sess}{i})
                    matfile = fullfile(mfdir,DSGN.pmods{sess}{i}{j});
                    matfile = [regexprep(matfile,'.mat$','') '.mat'];
                    if ~exist(matfile,'file'), error('%sNo such pmods file: %s',z,matfile); end
                    pmodinfo = load(matfile);
                    
                    for k = 1:numel(pmodinfo.name)
                        l = size(currpmods,2) + 1;
                        % initialize
                        currpmods(l) = newpmod;
                        currpmods(l).name = pmodinfo.name{k};
                        currpmods(l).param = pmodinfo.param{k};
                        currpmods(l).poly = pmodinfo.poly{k};
                    end
                end
            catch %#ok
                currpmods = emptypmod;
            end
        else
            currpmods = emptypmod;
        end
        
        for m = 1:numel(matfiles)
            condinfo = load(matfiles{m});
            
            % loop through all conditions in mat
            for j = 1:numel(condinfo.name)
                if ~numel(condinfo.onset{j})
                    names{session}{c} = condinfo.name{j}; %#ok
                    onsets{session}{c} = []; %#ok
                    durations{session}{c} = []; %#ok
                    pmods{session}{c} = emptypmod; %#ok
                else
                    % load name
                    names{session}{c} = condinfo.name{j}; %#ok
                    
                    % load onsets
                    if size(condinfo.onset{j},1) == 1 && size(condinfo.onset{j},2) > 1
                        condinfo.onset{j} = condinfo.onset{j}';
                    elseif size(condinfo.onset{j},1) > 1 && size(condinfo.onset{j},2) > 1
                        error('onsets field in %s stored as matrix (must be scalar or vector)',matfiles{m});
                    end
                    onsets{session}{c} = condinfo.onset{j}; %#ok
                    
                    % load durations
                    if size(condinfo.duration{j},1) == 1 && size(condinfo.duration{j},2) > 1
                        condinfo.duration{j} = condinfo.duration{j}';
                    elseif size(condinfo.duration{j},1) > 1 && size(condinfo.duration{j},2) > 1
                        error('durations field in %s stored as matrix (must be scalar or vector)',matfiles{m});
                    end
                    durations{session}{c} = condinfo.duration{j}; %#ok
                    
                    % add currently indicated pmods
                    pmods{session}{c} = currpmods; %#ok
                    
                    % load pmods from condition file                   
                    if isfield(condinfo,'pmod') && ~isempty(condinfo.pmod)
                        for k = 1:numel(condinfo.pmod.name)
                            l = size(pmods{session}{c},2) + 1;
                            pmods{session}{c}(l) = newpmod; %#ok
                            pmods{session}{c}(l).name = condinfo.pmod.name{k}; %#ok
                            pmods{session}{c}(l).param = condinfo.pmod.param{k}; %#ok
                            pmods{session}{c}(l).poly = condinfo.pmod.poly{k}; %#ok
                        end
                    end
                end
                c=c+1;
            end
        end
        
    end        
    
    % retrieve multiple regressors file
    if ~isempty(DSGN.multireg)
        multiregfile = fullfile(mfdir, DSGN.multireg);
%         if ~exist(multiregfile,'file'), error('> No such multiple regressors file : %s',multiregfile); end
        multipleregressors{session} = multiregfile; %#ok                
    else
        multipleregressors{session} = {}; %#ok
    end    
end

end


%%
function [runs3d names onsets durations pmods multipleregressors] = concatdata(DSGN,submodeldir,oldruns,oldruns3d,oldnames,oldonsets,olddurations,oldpmods,oldmultipleregressors,z)


emptyruns = find(cellfun('isempty',oldruns));
if DSGN.allowmissingfunc && ~isempty(emptyruns)        
    i=1;
    for c = 1:numel(DSGN.concatenation)
        concat{i} = []; %#ok
        for r = 1:numel(DSGN.concatenation{c})
            if ~any(DSGN.concatenation{c}(r) == emptyruns)
                concat{i} = [concat{i} DSGN.concatenation{c}(r)]; %#ok                
            end
        end
        if (numel(concat{i}) > 1), i=i+1; end            
    end
else
    concat = DSGN.concatenation;
end

for i = 1:numel(concat)
    fprintf('%s\tsession %d is run(s):',z,i);
    for j = 1:numel(concat{i})
        fprintf('  %3d', concat{i}(j));
    end
    fprintf('\n');
end

% initialize
runs3d = {};
names = {};
onsets = {};
durations = {};
pmods = {};
multipleregressors = '';

for sess = 1:numel(concat)  
    oldsess1 = concat{sess}(1);
    
    % concatenate functional data
    runs3d{sess} = []; %#ok
    for r = 1:numel(concat{sess})
        oldsess = concat{sess}(r);
        runs3d{sess} = [runs3d{sess}; cellstr(expand_4d_filenames(oldruns{oldsess}))]; %#ok
    end   
    
    starttime = 0;
    for r = 1:numel(concat{sess})
        oldsess = concat{sess}(r);
        for cond = 1:numel(oldonsets{oldsess})
            % names
            if r==1
                names{sess}{cond} = oldnames{oldsess}{cond}; %#ok
            elseif ~strcmp(names{sess}{cond},oldnames{oldsess}{cond})
                error(['Inconsistent names for condition ' num2str(cond) ' across sessions (e.g., ' names{sess}{cond} ', ' oldnames{oldsess}{cond} ')'])
            end
                
            % onsets
            if r==1, onsets{sess}{cond} = []; end %#ok
            onsets{sess}{cond} = [onsets{sess}{cond}; oldonsets{oldsess}{cond} + starttime]; %#ok
            
            % durations
            if r==1, durations{sess}{cond} = []; end %#ok
            if numel(olddurations{oldsess}{cond}) == 1
                % extend single duration across all onsets (in case single duration is different across runs being concatenated)
                olddurations{oldsess}{cond} = repmat(olddurations{oldsess}{cond},size(oldonsets{oldsess}{cond},1),1);
            end
            durations{sess}{cond} = [durations{sess}{cond}; olddurations{oldsess}{cond}]; %#ok            
            
            % pmods
            if ~numel(oldpmods{oldsess}{cond})
                pmods{sess}{cond} = oldpmods{oldsess}{cond}; %#ok
            else
                for p = 1:size(oldpmods{oldsess}{cond},2)
                    if r==1
                        pmods{sess}{cond}(p).poly = oldpmods{oldsess1}{cond}(p).poly; %#ok
                        pmods{sess}{cond}(p).name = oldpmods{oldsess1}{cond}(p).name; %#ok
                        pmods{sess}{cond}(p).param = []; %#ok
                    end
                    pmods{sess}{cond}(p).param = [pmods{sess}{cond}(p).param; oldpmods{oldsess}{cond}(p).param]; %#ok
                end
            end
        end
        starttime = starttime + (DSGN.tr * numel(oldruns3d{oldsess}));
    end
    
    % concatenate regressors
    multipleregressors{sess} = fullfile(submodeldir,sprintf('multireg_%d.mat',sess));
    newR = [];
    cri = {};
    for r = 1:numel(concat{sess})
        oldsess = concat{sess}(r);
        if ~isempty(oldmultipleregressors{oldsess})
            load(oldmultipleregressors{oldsess});
            oldR = R;            
        else
            oldR=[];
        end
        if ~isfield(DSGN,'customrunintercepts')
            % add intercept (ignore first one)
            if r>1
                oldR(:,end+1) = 1; %#ok
            end
        else
            % initialize
            if isempty(oldR)
                tmpn = nifti(oldruns{oldsess});
                cri{r} = zeros(size(tmpn.dat,4),numel(DSGN.customrunintercepts)); %#ok
            else
                cri{r} = zeros(size(oldR,1),numel(DSGN.customrunintercepts)); %#ok
            end
            
            for i = 1:numel(DSGN.customrunintercepts)
                if any(oldsess == DSGN.customrunintercepts{i})
                    cri{r}(:,i) = 1; %#ok
                end
            end
        end
        % add linear trend
        oldR(:,end+1) = scale([1:size(oldR,1)]'); %#ok
        
        % append to growing block diagonal nuisance matrix
        newR = blkdiag(newR,oldR);        
    end    
    R = [newR vertcat(cri{:})];
    save(multipleregressors{sess}, 'R');
end

% add rest of stuff
catruns = cell2mat(concat);
for r = 1:numel(oldruns)
    if ~any(catruns == r) && ~isempty(oldruns{r})
        fprintf ('%s\tsession %d is run(s): %3d\n',z,numel(runs3d)+1,r)        
        runs3d{end+1} = oldruns3d{r}; %#ok
        names{end+1} = oldnames{r}; %#ok
        onsets{end+1} = oldonsets{r}; %#ok
        durations{end+1} = olddurations{r}; %#ok
        pmods{end+1} = oldpmods{r}; %#ok
        multipleregressors{end+1} = oldmultipleregressors{r}; %#ok
    end
end

end


%%
function [names onsets durations pmods] = convert_to_single_trials(DSGN,oldnames,oldonsets,olddurations,oldpmods)

newpmod = struct('name', [], 'param', [], 'poly', []);
emptypmod = newpmod([]);

if isfield(DSGN,'singletrials') && numel(DSGN.singletrials) == 1
    for s = 2:numel(oldnames)
        DSGN.singletrials{s} = DSGN.singletrials{1};
    end
end

for s = 1:numel(oldonsets)
    o = 1;
    for c = 1:numel(oldonsets{s})
        if numel(olddurations{s}{c})==1
            olddurations{s}{c} = repmat(olddurations{s}{c},numel(oldonsets{s}{c}),1);
        end
        % note: look for better way to deal with sparse cell array!
        try 
            thiscond = logical(DSGN.singletrials{s}{c});
            if isempty(DSGN.singletrials{s}{c}), DSGN.singletrials{s}{c} = false; end
        catch %#ok
            thiscond = false;
        end
         
        if DSGN.singletrialsall || thiscond
            if sum(size(oldpmods{s}{c})==[0 0])~=2                                
                error('Sorry, single trials option is not currently compatible with pmods')
            end
            for t = 1:numel(oldonsets{s}{c})
                names{s}{o} = sprintf('%s_trial%04d',oldnames{s}{c},t); %#ok
                onsets{s}{o} = oldonsets{s}{c}(t); %#ok
                durations{s}{o} = olddurations{s}{c}(t); %#ok                
                pmods{s}{o} = emptypmod; %#ok
                o=o+1;
            end
        else
            names{s}{o} = oldnames{s}{c}; %#ok
            onsets{s}{o} = oldonsets{s}{c}; %#ok
            durations{s}{o} = olddurations{s}{c}; %#ok
            pmods{s}{o} = oldpmods{s}{c}; %#ok
            o=o+1;
        end
    end
end

end


%%
function [nonemptyruns nonemptyruns3d flatnames flatonsets flatdurations flatpmods conditions_by_run OPTS] = prep_for_canlab_spm_fmri_model_job(DSGN,OPTS,runs,runs3d,names,onsets,durations,pmods,z)

nonemptyruns = {};
nonemptyruns3d = {};
for i=1:numel(runs3d)
    if ~isempty(runs3d{i})
        nonemptyruns{end+1} = runs{i}; %#ok
        nonemptyruns3d{end+1} = runs3d{i}; %#ok
    end
end

flatnames = {};
flatonsets = {};
flatdurations = {};
flatpmods = {};
conditions_by_run = [];
newsess = 0;
for session = find(~cellfun('isempty',runs3d)) %1:numel(names)
    newsess = newsess+1;
    i=0;
    for cond = 1:numel(names{session})
        if isempty(onsets{session}{cond})
            if DSGN.allowemptycond
                fprintf('%sWARNING: no onsets in session %d, condition %d: %s\n',z,session,cond,names{session}{cond})
                continue
            else
                error('no onsets in session %d, condition %d: %s\n\t(to allow, set DSGN.allowemptycond)',session,cond,names{session}{cond})
            end
        end
        flatnames{end+1} = [names{session}{cond} ' ']; %#ok % space added to separate out tmods, pmods, basis functions, etc 
        flatonsets{end+1} = onsets{session}{cond}; %#ok
        flatdurations{end+1} = durations{session}{cond}; %#ok
        flatpmods{end+1} = pmods{session}{cond}; %#ok
        i=i+1;
    end
    conditions_by_run(newsess) = i; %#ok
end

switch DSGN.convolution.type
    case 'hrf'
        OPTS.modeljob = [OPTS.modeljob ',' '''hrf''' ',' num2str(DSGN.convolution.time) ',' num2str(DSGN.convolution.dispersion)];
    case 'fir'
        OPTS.modeljob = [OPTS.modeljob ',' '''fir''' ',' num2str(DSGN.convolution.windowlength) ',' num2str(DSGN.convolution.order)];
        if ~isfield(DSGN.convolution,'keepdurations') || ~DSGN.convolution.keepdurations
            % zero-out durations
            for c = 1:numel(flatdurations)
                flatdurations{c} = 0; %#ok
            end
        end
    otherwise
        error('Unrecognized convolution type: %s',DSGN.convolution.type)
end

if DSGN.ar1
    OPTS.modeljob = [OPTS.modeljob ',''AR(1)'''];
end

end


%%
function link_to_stat_maps(targdir)

startingdir = pwd;

cd(targdir)
load('SPM.mat')
% load('contrastnames.mat')

renamedir = fullfile(pwd, 'named_statmaps');
if exist(renamedir,'dir'), rmdir(renamedir,'s'); end
mkdir(renamedir);

for i=1:numel(SPM.xCon)
    mapname = SPM.xCon(i).name;
    % clean up the name for filename friendliness
    mapname = regexprep(mapname,'[()]','');
    mapname = regexprep(mapname,'[^0-9A-Za-z-_]','_');
    for stat = {'spmT' 'con'}
        for ext = {'img' 'hdr'}
            imgname = fullfile(pwd, sprintf('%s_%04d.%s',stat{1},i,ext{1}));
            linkname = fullfile(renamedir, [stat{1} '_' mapname '.' ext{1}]);
            eval(['!ln -v -s ' imgname ' ' linkname]);
        end
    end
end

cd(startingdir)

end
