function [infos] = canlab_glm_getinfo(modeldir,varargin)
% :SUBJECT LEVEL input:
% ::
%
%    INFO = canlab_glm_getinfo(spm_subject_level_directory, option, [n])
%
% Get information out of an spm subject level analysis.
%
% :Options: (each option can be called by the listed letter or word + a number, when noted)
%
%   **i' 'input':**
%        number of volumes and (first) volume name for each run
%
%   **'b' 'betas' [n]:**
%        beta names (for nth session)     
%
%   **'B' 'taskbetas' [n]:**
%        beta names (that didn't come from multiple regressors)  (for nth session)
%
%   **'c' 'cons' [n]:**
%        contrast names (for nth contrast)
%
%   **'C' 'conw' [n]:**
%        beta names and weights for contrasts (for nth con)
%
%   **'v' 'image' [n]:**
%        create figure of design matrix (for nth session)
%                            (design matrix is multiplied by 100 for visibility)
%                            (works well for multiple runs)
%
%   **'V' 'taskimage' [n]:**
%        same as 'image', but only for task betas
%
%   **'imagesc':**
%        same as 'image', but uses imagesc
%                            (works well for single runs)
%
%
% :GROUP LEVEL input:
% ::
%
%    INFO = canlab_glm_getinfo(robfit_group_level_directory, option, [n])
%
% Get information out of a robfit group level analysis.
%
% :Options: (each option can be called by the listed word or letter + a number, when noted)
%
%   Any of the subject level options can be used on a group level robfit
%   analysis by prefixing '1i' (output is generated based on the first input
%   to the first robust analysis).
%
%   Ex:
%   ::
%
%      canlab_glm_getinfo('second_level/model3','1iconw')
%
%
%   **'i' 'input' [n]:**
%        input contrasts by number and name (for nth analysis)
%
%   **'I' 'allinput' [n]:**
%        input images (for nth analysis)
%
%   **'m' 'model':**
%        weights by subject (i.e., directory containing input contrast images)
%
%   **'M' 'allmodels' [n]:**
%        weights and input images (for nth analysis)
%
% :Assumptions:
%    In some options, the first contrasts and group level analysis
%    directories are assumed to represent the rest, which may not be the
%    case.
%
% :Note:
%    group level options do not yet return a usable INFO struct.

if nargin ~= 2 && nargin ~= 3
    error('USAGE: INFO = canlab_glm_getinfo(spm_dir|robfit_dir,option,[n])')
end

opt = varargin{1};
if nargin == 3, s = varargin{2}; else s = 0; end

if ~exist(modeldir,'dir')
    error('No such directory: %s',modeldir)
else
    spmmat = fullfile(modeldir,'SPM.mat');
    setupmats = filenames(fullfile(modeldir,'robust[0-9][0-9][0-9][0-9]','SETUP.mat'));
    if exist(spmmat,'file')
        infos = get_subject_level_info(modeldir,opt,s);
    elseif ~isempty(setupmats)
        if regexp(opt,'^1i')
            load(setupmats{1});
            firstinput = fileparts(deblank(SETUP.files(1,:)));
            infos = get_subject_level_info(firstinput,regexprep(opt,'^1i',''),s);
        else
            infos = get_group_level_info(modeldir,opt,s);
        end            
    else
        error('%s is neither a subject level SPM directory nor a group level robfit directory',modeldir);
    end
end


end


function [infos] = get_subject_level_info(modeldir,opt,s)

load(fullfile(modeldir,'SPM'));

switch opt
    case {'c' 'con' 'cons'}
        if s
            infos.contrast_name{1} = SPM.xCon(s).name;
        else            
            infos.contrast_name = cellstr(strvcat(SPM.xCon.name));  %#ok
        end
        
        fprintf('CONTRASTS:\n')
        for i=1:numel(infos.contrast_name)
            fprintf('%4d %s\n',i,infos.contrast_name{i});
        end
        
    case {'C' 'conw' 'consw'}
        if s
            infos.contrast_name{1} = SPM.xCon(s).name;
            infos.beta_numbers{1} = find(SPM.xCon(s).c);
            infos.beta_names{1} = cellstr(strvcat(SPM.xX.name{infos.beta_numbers{1}})); %#ok
            infos.beta_weights{1} = SPM.xCon(s).c(infos.beta_numbers{1});
        else            
            infos.contrast_name = cellstr(strvcat(SPM.xCon.name)); %#ok            
            infos.contrast_num = 1:numel(infos.contrast_name);
            for i = 1:numel(infos.contrast_name);
                infos.beta_numbers{i} = find(SPM.xCon(i).c);
                infos.beta_names{i} = cellstr(strvcat(SPM.xX.name{infos.beta_numbers{i}})); %#ok
                infos.beta_weights{i} = SPM.xCon(i).c(infos.beta_numbers{i});
            end
        end
        
        for s = 1:numel(infos.contrast_name)
            fprintf('%3d: %s:\n',s,infos.contrast_name{s})
            for i=1:numel(infos.beta_names{s})
                fprintf('\t%10.5f\t%s\n',infos.beta_weights{s}(i),infos.beta_names{s}{i});
            end
        end
        
        
    case {'b' 'betas'}
        if s
            infos.beta_number = find(~cellfun('isempty',regexp(SPM.xX.name,sprintf('^Sn\\(%d\\)',s))))';
            infos.beta_name = cellstr(strvcat(SPM.xX.name{infos.beta_number})); %#ok
        else
            infos.beta_name = cellstr(strvcat(SPM.xX.name)); %#ok
            infos.beta_number = [1:numel(infos.beta_name)]'; %#ok
        end
        
        for i = 1:numel(infos.beta_name)
            fprintf('%4d %s\n',infos.beta_number(i),infos.beta_name{i});
        end
        
        
    case {'B' 'taskbetas'}
        n = cellfun('isempty',regexp(SPM.xX.name,' R[0-9]*$'))';
        if s
            n = n & ~cellfun('isempty',regexp(SPM.xX.name,sprintf('^Sn\\(%d\\)',s)))';
        end
        n = n & cellfun('isempty',regexp(SPM.xX.name,'constant$'))';
        infos.beta_number = find(n);
        infos.beta_name = cellstr(strvcat(SPM.xX.name{infos.beta_number})); %#ok
        
        for i = 1:numel(infos.beta_number)
            fprintf('%4d %s\n',infos.beta_number(i),infos.beta_name{i});
        end        
        
        
    case 'nsess'
        infos.nsess = numel(SPM.Sess);
        fprintf('%d\n',infos.nsess);
        
        
    case {'i' 'input' 'sess'}
        infos.nsess = numel(SPM.nscan);
        infos.nscan = SPM.nscan;
        
        n=1; for i=1:infos.nsess, first(i)=n; n=n+infos.nscan(i); end %#ok
        
        fprintf('%4s %5s  %s\n','sess','nscan','first_volume_name')
        for i = 1:numel(first)
            infos.first_volume_name{i} = SPM.xY.VY(first(i)).fname;
            fprintf('%4d %5d  %s\n',i,infos.nscan(i),infos.first_volume_name{i})
        end
        
        
    case {'v' 'image'}        
        if s
            canlab_glm_getinfo(modeldir,'betas',s);
            infos.X = SPM.xX.X(SPM.Sess(s).row,SPM.Sess(s).col);
        else
            canlab_glm_getinfo(modeldir,'betas');
            infos.X = SPM.xX.X;            
        end
        
        figure; image(infos.X * 40);
        colormap('Bone');
        
        
    case {'V' 'taskimage'}        
%         n = cellfun('isempty',regexp(SPM.xX.name,' R[0-9]*$'));
        if s
            n = canlab_glm_getinfo(modeldir,'taskbetas',s);
            %n = n .* ~cellfun('isempty',regexp(SPM.xX.name,['^Sn\(' num2str(s) '\)']));
            infos.X = SPM.xX.X(SPM.Sess(s).row,n.beta_number);            
        else
            n = canlab_glm_getinfo(modeldir,'taskbetas');            
            infos.X = SPM.xX.X(:,n.beta_number);            
        end        
        figure; image(infos.X * 400)
        colormap('Bone');
        
    case {'vsc' 'imagesc'}        
        spmlowerinfo(modeldir,'betas')
        if s
            canlab_glm_getinfo(modeldir,'betas',s)
            figure; imagesc(SPM.xX.X(SPM.Sess(s).row,SPM.Sess(s).col))
            colormap('Bone');
        else
            canlab_glm_getinfo(modeldir,'betas')
            figure; imagesc(SPM.xX.X)
            colormap('Bone');
        end
        
        
    otherwise
        error('Unrecognized option for spm lower level analyses: %s',opt)
end 

end

function [infos] = get_group_level_info(modeldir,opt,s)

infos = [];

% load SETUP structs
f = filenames(fullfile(modeldir,'robust[0-9][0-9][0-9][0-9]'));
for i=1:numel(f),
    load(fullfile(f{i},'SETUP'));
    SETUPS(i) = SETUP; %#ok
    clear SETUP
end

switch opt
    case {'i' 'input'}
        % load lower levels        
        evalc('tmp = canlab_glm_getinfo(fileparts(SETUPS(i).V(1).fname),''cons'');');
        connames = tmp.contrast_name;
        if ~s, s = 1:numel(SETUPS); end
        fprintf('%-12s %-12s %-s\n','dir','input_con','con_name')        
        for i = s          
            [ignore finf] = fileparts(SETUPS(i).V(1).fname); %#ok
            
            fprintf('%-12s',SETUPS(i).dir)
            fprintf(' %-12s',finf)
            if regexp(finf,'^con_')
                fprintf(' "%-s"',connames{str2num(regexprep(finf,'^con_',''))}) %#ok
            end
            fprintf('\n')
        end
    case {'I' 'allinput'}
        if ~s, s = 1:numel(SETUPS); end
        for i = s
            fprintf('%s:\n',SETUPS(i).dir)
            disp(SETUPS(i).files)
        end
    case {'m' 'model'}
        fprintf('weight input\n')
        for i = 1:numel(SETUPS(1).X)
            fprintf('%-7.3f %s\n',SETUPS(1).X(i),regexprep(SETUPS(1).V(i).fname,'/[^/]*$',''))
        end
    case {'M' 'allmodels'}
        if ~s, s = 1:numel(SETUPS); end
        for m = s
            fprintf('%s\n',SETUPS(m).dir)
            fprintf('weight input\n')
            for i = 1:numel(SETUPS(m).X)
                fprintf('%-7.3f %s\n',SETUPS(m).X(i),SETUPS(m).V(i).fname)
            end
        end
        
    otherwise
        error('Unrecognized option for robfit analyses: %s',opt)
end



end
        
        
