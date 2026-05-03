% canlab_glm_publish_group_levels Publish group-level robust regression results via robust_results_batch.
%
% :Usage:
% ::
%
%     % Run from a directory containing robust* subdirectories
%     canlab_glm_publish_group_levels
%
% Child script of canlab_glm_publish. Runs robust_results_batch in a
% publishable manner so the results can be rendered as an HTML report
% via MATLAB's publish().
%
% Behavior:
%
%   - If a variable EXPT exists in the workspace, contrast names
%     and robust#### directories are taken from EXPT.SNPM.
%   - Otherwise the script discovers robust[0-9][0-9][0-9][0-9]
%     directories in the current working directory; if none are found,
%     the current directory itself is treated as a single robust*
%     directory. Contrast names are read from each directory's
%     SETUP.mat / SPM.mat.
%   - A small wrapper script robfit_results.m is written and then
%     published. The output HTML report is placed in
%     Robust_Regression_Results_html under the working directory.
%
% Workspace variables thresh and size (passed to
% robust_results_batch) default to [.001 .005 .05] and [5 1 1]
% respectively if not already defined.
%
% :Inputs:
%
%   None. This file is a script and operates on the current working
%   directory and the workspace variables EXPT, thresh, and
%   size if they exist.
%
% :Outputs:
%
%   No MATLAB outputs. Side effects:
%
%   - Creates Robust_Regression_Results_html directory with the
%     published HTML report.
%   - Creates and then deletes the helper script robfit_results.m
%     in the working directory.
%
% :See also:
%   - canlab_glm_publish_subject_levels
%   - canlab_glm_group_levels
%   - robust_results_batch
%   - publish

% child script of canlab_glm_publish
% runs robust_results_batch in publishable manner
% runs in directory containing robust* directories

%% prep
z = '________________________________________________________________________________________________________';
assignin('base','z',z);
if ~exist('thresh','var'), thresh = [.001 .005 .05]; assignin('base','thresh',thresh); end
if ~exist('size','var'), size = [5 1 1]; assignin('base','size',size); end

if exist('EXPT','var')
    for i = 1:length(EXPT.SNPM.P)
        robregdirs{i} = fullfile(pwd,sprintf('robust%04d',i));
        connames{i} = EXPT.SNPM.connames(i,:);
    end    
else
    d = filenames('robust[0-9][0-9][0-9][0-9]','absolute');
    if isempty(d), d{1} = pwd; end
    for i = 1:numel(d)
        robregdirs{i} = d{i};
        load(fullfile(robregdirs{i},'SETUP.mat'));
        load(fullfile(fileparts(SETUP.files(1,:)),'SPM.mat'));
        c = regexprep(SETUP.files(1,:),'.*con_0*([0-9]*)\.nii','$1');
        connames{i} = SPM.xCon(str2num(c)).name;
    end
end


%% write script
scriptfile = 'robfit_results.m';
scriptfileabs = fullfile(pwd,scriptfile);
fid = fopen(scriptfile,'w');

for i = 1:length(robregdirs)        
    fprintf(fid,'%%%% contrast(%d) %s\n',i,connames{i});
    fprintf(fid,'cd(''%s'')\n',robregdirs{i});
    fprintf(fid,'fprintf(''%%s\\n%%s\\n%%s\\n%%s\\n%%s\\n%%s\\n'',z,z,''%s'',''%s'',z,z)\n',connames{i},robregdirs{i});
    
    fprintf(fid,'try\n');
    maskimg = fullfile(robregdirs{i}, 'rob_tmap_0001.img'); % any image in space with non-zero vals for all vox would do
    %fprintf(fid,'robust_results_batch(''thresh'', thresh, ''size'', size, ''prune'', ''mask'', ''%s'');\n',maskimg);
    fprintf(fid,'robust_results_batch;\n');
    
    fprintf(fid,'catch exc\n');
    fprintf(fid,'disp(getReport(exc,''extended''))\n'); 
    fprintf(fid,'end\n');
end
fprintf(fid,'close all\n');

fclose(fid);


%% publish script
outputdir = fullfile(pwd, 'Robust_Regression_Results_html');
mkdir(outputdir);

p = struct('useNewFigure', false, 'maxHeight', 1500, 'maxWidth', 1200, ...
    'outputDir', outputdir, 'showCode', false);

fout = publish(scriptfile, p);
fprintf('Created robfit results directory:\n%s\nHTML report: %s\n', outputdir, fout);


%% clean up
close all
delete(scriptfileabs)
