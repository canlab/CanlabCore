function [robfitstatus] = canlab_glm_group_levels_run1input(wd, c)
% child process of canlab_glm_group_levels
% (see canlab_glm_README.txt for an overview)
%
% ..
%    Copyright (C) 2013  Luka Ruzic
% ..

load(fullfile(wd,sprintf('env_%04d',c)));
%% PREP

% diaryname = fullfile(wd,sprintf('diary_%04d.log',c));
cd(grpmodeldir)

robfitstatus = 0; %#ok
% grfstatus = 0;

%% robfit
cmd = 'robfit(EXPT, includedcons(c), 0, EXPT.mask)';
% diary(diaryname), fprintf('> %s\n',cmd); diary off
fprintf('> %s\n',cmd);
try    
    eval(cmd)
    robfitstatus = 1;
catch exc
    if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc);
    else fprintf('> %s\n',getReport(exc,'extended')); end
    %     else diary(diaryname), fprintf(getReport(exc,'extended')), diary off; end
    robfitstatus = -1;
end


%% grf  % having problems with this running in parallel (display problems?)
% if OPTS.run_grf
%     robdir = filenames(fullfile(grpmodeldir,sprintf('robust%04d',includedcons(c))),'char','absolute');
%     if ~isempty(robdir)
%         announce_string('SIGNIFICANT CLUSTER SIZE ESTIMATION using GRF')
%         try
%             load(fullfile(robdir,'SETUP.mat'));
%             
%             cmd = 'sigclext = estimate_cluster_extent(.05, pthresh, SETUP.files);';
%             fprintf('> %s\n',cmd)
%             eval(cmd);
%             close all
%             
%             fout = fullfile(robdir,'significant_cluster_extents_grf.txt');
%             dlmwrite(fout,[pthresh' sigclext(:,1)],'precision','%g','delimiter',' '); %#ok
%        
%         grfstatus = 1;
%         catch exc
%             if OPTS.nocatch, cd(STARTINGDIR); rethrow(exc)
%             else fprintf('> %s\n',getReport(exc,'extended')); end
%             grfstatus = -1;
%         end
%     end
% end

end



function announce_string(string)

s = sprintf('--  %s  --',string);
l = regexprep(s,'.','-');
fprintf('> \n> \n> \n> %s\n> %s\n> %s\n> \n',l,s,l);

end
