function canlab_glm_publish(varargin)
% :Usage:
% ::
%
%     canlab_glm_publish(directory_specifications [options])
%
% :Directory Specification:
%
%   **'s', dirs:**
%        Generates HTML reports of the output from scn_spm_design_check for
%        directories in cell array dires (string allowable for single dir).
%        If a directory is a subject level analysis, the HTML will be generated 
%        for that subject in the analysis directory.
%        Ex:
%        ::
%
%        canlab_glm_publish('s',{'model1/1011' 'model1/1014'})
%
%        If a directory contains subject level analyses, an HTML will be
%        generated with output for each subject level analysis.
%        Ex:
%        ::
%
%        canlab_glm_publish('s','model1')
%
%        ASSUMPTION: lower level analyses contain an SPM.mat file.
%
%
%   **'g', dirs:**
%        For each "robfit directory" in cell array dirs, will run robust_results_batch 
%        on all contrast directories (e.g., robust0001/) (string allowable for single dir).
%        ("robfit directories" contain robfit contrast directories (like robust0001))
%        EITHER directories contain EXPT.mat files (see help robfit)
%        OR an EXPT struct is loaded in the workspace and a single directory is specified
%        OR will do best to figure out info normally contained in EXPT.mat
%        Ex:
%        ::
%
%        canlab_glm_publish('g', {'group_n35' 'group_anxiety_n35' 'group_sadness_n35'})
%
% :Note: directory paths may be absolute or relative (to working directory)
%
%
% :Options:
%
%   **'t', {[pthresh clustersize] ...}:**
%        Use the paired voxelwise_pthresh and minimum_cluster_size thresholds 
%        with which to produce robfit results maps.
%        This option must follow immediately after a 'g' option (see above) and
%        will only apply to the analyses specified in that option.
%        ONLY applies to robfit directories (no bearing on lower level design checks)
%
%        DEFAULT: {[.001 5] [.005 1] [.05 1]}
%        Ex:
%        ::
%
%        canlab_glm_publish('g', pwd, 't', {[.001 1] [.005 10] [.05 10] [.01 25]})
%
%
%   **'email', address:**
%        send notification email to address when done running
%        Ex:
%        ::
%
%        canlab_glm_publish('g', pwd, 'email', 'ruzic@colorado.edu')
%


thresh = [.001 .005 .05];
size = [5 1 1];

%% setup
STARTINGDIR = pwd;

%% parse arguments
    
i=1;
while i<=numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'s' 'l'}
                i=i+1;
                if isstr(varargin{i})
                    varargin{i} = {varargin{i}};
                elseif ~iscellstr(varargin{i})
                    error('[Argument to ''s'' must a be string or cell array.]')
                end
                spmlower_report(STARTINGDIR,varargin{i});
                
            case {'g' 'r'}                
                i=i+1;
                if isstr(varargin{i})
                    varargin{i} = {varargin{i}};
                elseif ~iscellstr(varargin{i})
                    error('[Argument to ''g'' must a be string or cell array.]')
                end
                robfitdirs = varargin{i};
                
                if i < nargin && strcmp(varargin{i+1},'t')
                    i=i+2;
                    t = cell2mat(varargin{i});
                    thresh = t(1:2:end);
                    size = t(2:2:end);
                end
                
                robfit_report(STARTINGDIR,robfitdirs,thresh,size);       
                
            case {'email'}
                i=i+1;
                address = varargin{i};
                
            otherwise
                error(['Unrecognized argument: ' varargin{i}])
                
        end
    else
        error(['Unrecognized argument: ' varargin{i}])
    end
    i=i+1;
end

%% clean up
cd(STARTINGDIR)

%% email notification
if exist('address','var')
    try
        [ignore output] = system(sprintf('printf "canlab_glm_publish has finished running.\n" | mail -v -s "canlab_glm_publish done" %s',address)); %#ok
    catch  %#ok
        fprintf('The notification email failed to send.\n');
    end
end

end


function spmlower_report(STARTINGDIR,analysisdirs)

for i = 1:numel(analysisdirs)
    if regexp(analysisdirs{i},'^/')
        targdir = analysisdirs{i};
    else
        targdir = fullfile(STARTINGDIR, analysisdirs{i});
    end
    
    if ~exist(targdir,'dir')
        warning('No such directory: %s', targdir)
        continue
    end
    
    
    cd(targdir);
    
    outputdir = fullfile(targdir, 'design_review_html');
    mkdir(outputdir)
    
    p = struct('useNewFigure', false, 'maxHeight', 1500, 'maxWidth', 1200, ...
        'outputDir', outputdir, 'showCode', false);
    
    fout = publish('canlab_glm_publish_subject_levels.m', p);
    fprintf('Created subject level design review:\n\t%s\n', fout);
end

end



%% robfit group analyses (robust_results_batch)
function robfit_report(STARTINGDIR,robfitdirs,thresh,size)

assignin('base','thresh',thresh)
assignin('base','size',size)

for i = 1:numel(robfitdirs)
    if regexp(robfitdirs{i},'^/')
        targdir = robfitdirs{i};
    else
        targdir = fullfile(STARTINGDIR, robfitdirs{i});
    end
    
    if ~exist(targdir,'dir')
        fprintf('*** WARNING: skipping robfit directory (not found): %s\n', targdir)
        continue
    end
    
    cd(targdir);
    canlab_glm_publish_group_levels
end

end
