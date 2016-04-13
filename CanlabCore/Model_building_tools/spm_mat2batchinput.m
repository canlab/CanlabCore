function [imgs, TR, scanspersess, names, onsets, durations] = spm_mat2batchinput(SPM)
% Extract info from SPM structure and format in batch input mode
%
% :Usage:
% ::
%
%     [imgs, TR, scanspersess, names, onsets, durations] = spm_mat2batchinput(SPM)
%
% :Input into SPM8::
%
%   tmod = time modulator (number per cell, order, linear = 1)
%
% For parametric modulation include a structure array, which is up to 1 x n in size, called pmod. n must be less
% than  or  equal  to  the  number  of  cells in the names/onsets/durations cell arrays. The structure array pmod
% must  have the fields: name, param and poly.  Each of these fields is in turn a cell array to allow the inclusion
% R for regressors of no interest
%
% NOTES:
% This is a preliminary version of this function and does not extract
% covariates, parametric modulators, etc.
%
% Returns output in the format you would use to enter into the SPM GUI to
% set up a model. Same format is used by canlab_spm_fmri_model_job.m
%
% ..
%    Tor Wager, Feb 2011
%    Could be extended to include nuisance covs, parametric modulators, etc.
%    Right now this is very basic.
% ..

if iscell(SPM.xY.P) % check input image names
    imgs = char(SPM.xY.P{:});
else
    imgs = SPM.xY.P;
end

TR = SPM.xY.RT;
scanspersess = SPM.nscan;

[names, onsets, durations, covs] = deal({});

% list regressor names
for sess = 1:length(SPM.Sess)
    
    names = [names cat(1, SPM.Sess(sess).U(:).name)'];
    %disp(sprintf('Session %3.0f', sess));
    
end


% check onsets
% 0 should be the start of the first image acquisition

for sess = 1:length(SPM.Sess)
    
    %disp(sprintf('Session %3.0f', sess));
    
    for evt = 1:length(SPM.Sess(sess).U)
        %         disp(['Session ' num2str(sess) ' event type: ' SPM.Sess(sess).U(evt).name])
        %         [SPM.Sess(sess).U(evt).ons  SPM.Sess(sess).U(evt).durations]
        
        onsets = [onsets {SPM.Sess(sess).U(evt).ons}];
        durations = [durations {SPM.Sess(sess).U(evt).dur}];
        
        
    end
    
    
end


