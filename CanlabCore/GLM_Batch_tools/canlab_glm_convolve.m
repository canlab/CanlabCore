% X = canlab_glm_convolve(stimObj, dsgn, nscan, basis)
%
% This function takes an spm style conditions object (see below), a canlab
% DSGN object (see below), a frame count and a basis function type and
% returns a timeseries vector corresponding to onsets and durations
% specified in stimObj and any modulations specified in dsgn, convolved
% with the basis function specified by basis. Useful if you want to treat
% experimental factors as confounds in a multiple regressors matrix instead
% of treating them as events in SPM. Primarly desinged to allow for
% multiple different basis functions to be used for different design
% factors (e.g. spline interpolation for factors of interest and canonical
% interpolation for confounding factors; here you would implement the
% confounding factors as multiple regressor columns instead after
% convolution).
%
% Input ::
%
%   stimObj     - an SPM-stype condition object. see
%                   canlab_glm_subject_levels('dsgninfo') under 
%                   DSGN.conditions for an example
%
%   dsgn        - see canlab_glm_subject_levels('dsgninfo') for how to
%                   format this object.
%
%   nscan       - number of TRs in the scan for which stimObj applies.
%
%   basis       - currently only 'hrf' is supported, although others could
%                   be implemented easily enough. See spm_get_bf.m for some
%                   example implementations of other basis functions. The
%                   implementation here follows that function closely.
%
% Output ::
%
%   X           - A timeseries vector corresponding to the stimuli
%                   specified in stimObj.
%
% Written by Bogdan Petre, 5/19/2021

function X = canlab_glm_convolve(stimObj, dsgn, nscan, basis)
% we're going to manually convolve stimFactor 
spm('defaults','fmri')
spm_jobman('initcfg')

fMRI_T = spm_get_defaults('stats.fmri.t');
if ismember('fmri_t',fieldnames(dsgn))
	if ~isempty(dsgn.fmri_t), fMRI_T = dsgn.fmri_t; end
end

T = dsgn.tr;
dt = dsgn.tr/fMRI_T;

switch basis
    case 'hrf'
            [bf, p] = spm_hrf(dt,[],fMRI_T);

        % time derivative (untested)
        try 
            if logical(dsgn.convolution.derivs(1))
                dp       = 1;
                p(6)     = p(6) + dp;
                D        = (bf(:,1) - spm_hrf(dt,p,fMRI_T))/dp;
                bf       = [bf D(:)];
                p(6)     = p(6) - dp;
            end
        end

        % dispersion derivative (untested)
        try 
            if logical(dsgn.convolution.derivs(2))
                dp   = 0.01;
                p(3) = p(3) + dp;
                D    = (bf(:,1) - spm_hrf(dt,p,fMRI_T))/dp;
                bf   = [bf D(:)];
            end
        end
    otherwise
        error('Only Canonical HRF convolution is currently supported');
end

bf = spm_orth(bf);

U = [struct(...
    'name','',...
    'ons',[],...
    'dur',[],...
    'orth', 1,...
    'P', struct('name','none','h',0))];


% rating and cue convolution
U(1).name = {'manualConvFactor'};
U(1).ons = stimObj.onset{1}(:);
U(1).dur = stimObj.duration{1}(:);

% this borrws from spm_fMRI_design.m which implements the
% convolution when spm is invoked directly.
xBF = struct('T',fMRI_T,'dt',dt,'UNITS','secs');
this = struct('Sess',struct('U',U), 'nscan', nscan, 'xBF', xBF);
U = spm_get_ons(this,1);
[X,~,~] = spm_Volterra(U, bf, 1);

if ~isempty(X)
    X = X((0:(this.nscan - 1))*fMRI_T + dsgn.fmri_t0 + 32,:);
end

end