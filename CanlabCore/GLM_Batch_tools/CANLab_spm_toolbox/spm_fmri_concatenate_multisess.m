function spm_fmri_concatenate_multisess(P, scans)
% Adjust an SPM.mat which has concatenated sessions.
% FORMAT spm_post_concatenate(P, scans)
% Session regressors are added and the high-pass filter and non-sphericity
% estimates adjusted as if sessions are separate.
%
% Unlike spm_fmri_concatenate_multisess this script should handle multisession
% designs with mixed concatenation/block diagonal designs. It has been
% tested for 4 scans of 1200 TRs with the first 2 and last 2 concatenated.
% YMMV.
% Bogdan Petre
% 11/17/2024
%
% P     - filename of the SPM.mat file to adjust
% scans - [1 x n] vector with the original number of scans in each session
%
% The expected workflow is:
%
% 1. Manually specify a GLM with timeseries and onsets concatenated
% 2. Run spm_post_concatenate on the saved SPM.mat.
% 3. Estimate the SPM.mat in the normal way.
%
% Tips:
%
% - The BOLD-response may overhang from one session to the next. To reduce
%   this, acquire additional volumes at the end of each session and / or
%   add regressors to model the trials at the session borders.
%__________________________________________________________________________
% Copyright (C) 2015-2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin & Peter Zeidman
% $Id: spm_fmri_concatenate.m 7018 2017-02-15 13:36:48Z guillaume $


%-Input parameters
%==========================================================================

%-Get SPM.mat and number of scans per session
%--------------------------------------------------------------------------
if nargin < 2
    error('Please specify number of scans per run/session')
end

% Load SPM.mat
%--------------------------------------------------------------------------
SPM = load(P);

% Validate
%--------------------------------------------------------------------------
try
    SPM = SPM.SPM;
catch
    error('Input file is not a valid SPM.mat.');
end

if sum(scans) ~= sum(SPM.nscan)
    error('Number of scans does not match.');
end

% Create backup
%--------------------------------------------------------------------------
copyfile(P, fullfile(spm_file(P,'fpath'),'SPM_backup.mat'));


%-Session-specific whitening and filtering
%==========================================================================

SPM.nscan = scans;

% Session-specific grand mean scaling (see spm_fmri_spm_ui.m)
%--------------------------------------------------------------------------

% Session effects (see spm_fMRI_design.m)
%--------------------------------------------------------------------------
n_sess = length(SPM.Sess);
if numel(SPM.xX.iB) == n_sess && all(SPM.xX.iB == (size(SPM.xX.X,2) - fliplr(0:n_sess-1)))
    Xb = [];
    Bn = {};
    for s=1:numel(SPM.nscan)
        Xb      = blkdiag(Xb,ones(SPM.nscan(s),1));
        Bn{s}   = sprintf('Sn(%i) constant',s);
    end
    n_intercepts = length(SPM.Sess);
    SPM.xX.X    = [SPM.xX.X(:,1:(end-n_intercepts)) Xb];
    SPM.xX.iB   = SPM.xX.iB:(SPM.xX.iB+size(Xb,2)-1);
    SPM.xX.name = {SPM.xX.name{1:end-n_intercepts} Bn{:}};

    % pop off any constants from session specific designs
    % nipype adds these automatically
    sess_ind = fliplr(1:length(SPM.Sess));
    for i = 1:length(sess_ind)
        s = sess_ind(i);
        sess_col = SPM.Sess(s).col;
        const_col = sess_col(contains(SPM.xX.name(sess_col),'constant') | ...
            ~cellfun(@isempty, regexp(SPM.xX.name(sess_col), 'UR\d+')));
        SPM.xX.X(:,const_col) = [];
        SPM.xX.name(const_col) = [];
        SPM.xX.iC(const_col) = [];
        while length(const_col) > 0
            SPM.xX.iC(SPM.xX.iC >= max(const_col)) = SPM.xX.iC(SPM.xX.iC >= max(const_col)) - 1;
            SPM.xX.iB(SPM.xX.iB >= max(const_col)) = SPM.xX.iB(SPM.xX.iB >= max(const_col)) - 1;
            
            % adjust indices of higher order sessions
            for j = s:length(SPM.Sess)
                this_sess_col = SPM.Sess(j).col;
                bad = ismember(this_sess_col, const_col);
                
                SPM.Sess(j).col(bad) = [];
                SPM.Sess(j).col(SPM.Sess(j).col > const_col) = SPM.Sess(j).col(SPM.Sess(j).col > const_col) - 1;

                if ~isempty(SPM.Sess(j).C.name), SPM.Sess(j).C.name(bad) = []; end
                if ~isempty(SPM.Sess(j).C.C), SPM.Sess(j).C.C(:,bad) = []; end
            end

            const_col(const_col == max(const_col)) = [];
        end
    end
end

% High-pass filter (see spm_spm.m)
%--------------------------------------------------------------------------
s = cumsum([0 SPM.nscan]);

for i=1:numel(SPM.nscan)
    K(i) = struct('HParam', SPM.xX.K(1).HParam,...
                  'row',    s(i) + (1:SPM.nscan(i)),...
                  'RT',     SPM.xY.RT);
end
SPM.xX.K = spm_filter(K);

% Temporal non-sphericity (see spm_spm.m)
%--------------------------------------------------------------------------
switch lower(SPM.xVi.form)
    case {'ar(1)','ar(0.2)'}
        SPM.xVi.Vi   = spm_Ce('ar',SPM.nscan,0.2);
        SPM.xVi.form = 'AR(0.2)';
    case 'fast'
        SPM.xVi.Vi   = spm_Ce('fast',SPM.nscan,SPM.xY.RT);
    case {'i.i.d', 'none'}
    otherwise
        warning('Unhandled temporal non-sphericity.');
end

% Set back nscan as if single session
%--------------------------------------------------------------------------
%SPM.nscan = sum(SPM.nscan);


%-Save
%==========================================================================
save(P, 'SPM', spm_get_defaults('mat.format'));
disp('SPM.mat adjusted for sessions: Please estimate to apply changes.');
