function delta = onsets2delta(ons, len)
% Builds high-res delta function, given cell array of onset times
%
% :Usage:
% ::
%
%     delta = onsets2delta(ons, [length: num rows in original units])
%
% Tor Wager, 2 / 24 / 04, update 10 / 12 / 10
%
% :Inputs:
%
%   **ons:**
%        onsets for each of a series of conditions
%
%        One cell per condition per session, e.g., ons{1} = [24 27 29 44]';
%
%        Units are arbitrary (e.g., TRs or seconds)
%
%        Onsets are assumed to start at time 0 (0 is start of run/session)
%
%        e.g., from an SPM.mat fmri design structure, for one session:
%        ::
%
%            ons = cell(1, nconds);
%            [ons{:}] = deal(reportmod_model.Sess(1).U(:).ons)
%
%   **len:**
%        optional: number of rows in original units.  Useful for making 
%      a design matrix with the right number of rows after convolution and downsampling
%
% :Output:
%
%   **delta:**
%        indicator matrix of vectors for each condition with 1/0 for each onset
%
%        Type is logical; you may want to do double(delta) before operating
%
%        Resolution of high-res delta functions = original units (secs or TRs) * 16
%
%        The number of rows is the max onset + 1, times 16
%
% For what to do with output, see:
% getPredictors : for design-matrix building
% downsample_canlab : for downsampling to TR/secs
%
% :See also: ONSETS2FMRIDESIGN and object-oriented fmri_model object
% (methods: build, etc.)
%
%
% ..
%    Revision: tor: 10/12/10, for integration with object-oriented fmri_model
% ..


% ..
%    Defaults
% ..

res = 16;   % resolution of high-res delta functions = orig. units * 16

if nargin == 1
    len = res * (max(cat(1, ons{:})) + 1);
else
    len = res * len;
end

% base indicator: zeros
cf = false(len,1);

nconds = length(ons);
delta = false(len, nconds); %cell(1, length(ons));

% ----------------------------------------------
% Error checking
% ----------------------------------------------
if ~iscell(ons), error('ons must be cell array'); end

if any(diff(cat(1, ons{:})) == 0)
    warning('onsets2delta: repeated onsets', 'Some onsets are identical');
end

% ----------------------------------------------
% Build
% ----------------------------------------------
for i = 1:nconds
    
    cf2 = cf;
    
    cf2(round(ons{i}*res) + 1) = true;   % first TR is time 0, element 1
    
    if length(cf2) > len, error('Some onsets occur after the end of the session. Mismatch in units??'); end
    
    delta(:, i) = cf2;
end

return

