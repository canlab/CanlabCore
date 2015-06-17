function [P] = hrf_power(TR,hrfsamples,delta,contrasts,varargin)
% [P] = hrf_power(TR,hrfsamples,delta,contrasts,[S],[X])
%
% Given a delta (indicator) matrix of onset times and contrasts, 
% this function generates an FIR deconvolution matrix and tests the 
% power of the design with pre-set assumptions.
%
%
% Tor Wager, 2 / 24 / 04
%
% Inputs:
% TR, rep time in s (sampling rate)
% hrfsamples, time points to estimate hemodynamic response
% delta, cell array containing onsets for each condition, sampled at TR
% contrasts, contrast matrix with one row per contrast 
%   ! first contrast is used to develop true response !
%   ! events with positive contrast weights are "active"
% S, optional high/low-pass filtering matrix
% X, optional pre-built FIR model matrix, for computational efficiency
%
% Outputs:
% P structure, with a number of fields:
% X2 = FIR deconvolution model
% ipow = individual subject 80% power, default reference effect
% gpow = group (default N) 80% power
% truehrf = hemodynamic true response for events of interest
%
% use with xpower.m and genMseq.m
%
% examples:
% delta is a cell array or matrix of 0's and 1's, two columns
% P = hrf_power(2,15,delta,[1 -1]);   est. response for 30 s

if isempty(hrfsamples), hrfsamples = round(30 ./ TR); end
P.hrfsamples = hrfsamples;
P.delta = delta;
P.contrasts = contrasts;

if length(varargin) > 0, P.S = varargin{1}; else P.S = []; end

% make design matrix
if length(varargin) > 1, 
   X2 = varargin{2}; 
else
    [X2] = tor_make_deconv_mtx3(delta,hrfsamples,1);
end

if ~isempty(P.S), X2 = P.S * X2; end
P.X2 = X2;

% Apply contrasts to HRF and build 'true' response

hrf2 = spm_hrf(TR); hrf2 = hrf2 ./ max(hrf2);
if hrfsamples > length(hrf2), hrf2 = [hrf2' zeros(1,hrfsamples-length(hrf2))]'; end

[P.hrfcontrasts,P.truehrf] = hrf_contrasts(contrasts,hrfsamples,hrf2);

% check
if size(P.hrfcontrasts, 2) > size(X2, 2)
    disp('Size mismatch!')
    P.hrfcontrasts = P.hrfcontrasts(:, 1:size(X2, 2));
    P.truehrf = P.truehrf(:, 1:size(X2, 2));
end

% calculate power
[P.ipow,P.gpow,P.PowerParams] = xpower(P.X2,P.hrfcontrasts,P.truehrf,[],[],[]);

return



function [hrfcon,truehrf] = hrf_contrasts(contrasts,hrfsamples,hrf)

hrfcon = []; truehrf = [];
if length(hrfsamples) < length(contrasts), hrfsamples = repmat(hrfsamples,1,length(contrasts));,end

for i = 1:length(contrasts)
    tmp = eye(hrfsamples(i)) .* contrasts(i);
    hrfcon = [hrfcon tmp];
    truehrf = [truehrf hrf(1:hrfsamples(i)) .* (contrasts(i) > 0)];
end
truehrf = truehrf(:);
if size(truehrf,1) > size(truehrf,2), truehrf = truehrf';,end

return



% extra code - old way

P.hrf = hrf2;

truehrf = []; conhrfmult = [];
tt = zeros(1,P.hrfsamples); 
for i = 1:size(P.delta,2),
   if P.contrasts(1,i)>0, 
        tt2=(tt+1).*hrf2(1:length(tt))';,
        conhrfm = tt + 1;
    else,
        tt2=tt;,
        conhrfm = tt; 
    end,
    truehrf=[truehrf tt2];,
    conhrfmult = [conhrfmult conhrfm];
end
P.hrfcontrasts = eye(size(X2,2) - 1);
P.hrfcontrasts(find(conhrfmult==0),:) = []; % contrasts only for HRFs of interest


