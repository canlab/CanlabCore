function [S,Vi,svi,KH] = getSmoothing(HPlength,LPsmooth,TR,numsamps,xc,varargin)
% [S,Vi,svi,KH] = getSmoothing(HPlength,LPsmooth,TR,numsamps,xc,[sessions])
%
% Inputs
% 	HPlength: 	high-pass filter length in s, or [] for none
%	LPsmooth: 	low-pass filter with hrf, 1 (yes) or 0 (no)
%	TR:			repetition time
%	numsamps:	number of scans in session
%	xc:			custom autocorrelation function, [] for white noise, 'auto'
%	            for canonical 1/f function
%
% Outputs
%	S			smoothing filter, apply using S * X
%	Vi			intrinsic autocorrelation matrix
%	svi			S * Vi * S';
%
% Tor Wager, 11/17/01, last modified 2/24/04 to speed processing with < 3
% output arguments.


% ----------------------------------------------------------------
% * get smoothing matrix and autocorrelation matrix
% ----------------------------------------------------------------
if isempty(HPlength),HPlength = 'none';,end
if isempty(LPsmooth),LPsmooth = 0;,end

% disp(['	...setting smoothing, hrf, and autocorrelation matrices, HPlength = ' num2str(HPlength) ', Smoothing = ' num2str(LPsmooth)])

%use spm_filter
if (strcmp(HPlength,'none') | isempty(HPlength)) & LPsmooth					% LP only
	%disp('		...LP only')
	[S,KL] = use_spm_filter(TR,numsamps,'hrf','none',[]);
elseif LPsmooth																% HP and LP

	[S,KL,KH] = use_spm_filter(TR,numsamps,'hrf','specify',HPlength);
	%disp('		...HP and LP')
elseif (strcmp(HPlength,'none') | isempty(HPlength)) & ~LPsmooth			% neither
	%disp('		...no filtering')
	%S = eye(numsamps);
	S = [];
else [S,KL,KH] = use_spm_filter(TR,numsamps,'none','specify',HPlength);		% HP only
	%disp('		HP only')
end

if nargout > 1

    if isempty(xc),
        Vi = eye(numsamps);,disp('Using white noise autocorrelation')
    elseif strcmp(xc,'auto')
        disp('Using canonical 1/f autocorrelation function')
        [xc,Vi] = canonical_autocorrelation(TR,numsamps);
    else 
        Vi = getv('make',xc,numsamps);
    end

    if nargout > 2
        
        if isempty(S)
	        svi = Vi;
        else
	        svi = S * Vi * S';
        end
    
    end
    
end

if nargout > 3 & length(varargin) > 0
% build model matrix for individual sessions
x = intercept_model(nsess);
end

    
    
end



