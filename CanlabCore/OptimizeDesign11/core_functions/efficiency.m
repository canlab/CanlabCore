function varargout = efficency(X,P)
% function [con_eff,hrf_eff,X,hrfX,P] = efficency(X,PARAMS)
% a useful summary function for a design, incorporating tests of contrasts and of HRF shape estimation
%
% X
%   either a design vector (condition function) with integers indicating onsets of conditions
%   or an already-constructed design matrix
%
% PARAMS structure contains inputs
% Optional fields:
%   contrasts   required for testing contrast efficiency; see calcEfficiency.m
%   Vi          intrinsic autocorrelation matrix; empty assumes identity matrix (independence)
%   S           smoothing matrix; empty for no smoothing
%   HRF         canonical hemodynamic response function; empty assumes SPM's canonical HRF
%               should be in resolution of .1 s per element
%   nonlint     saturation threshold (ceiling) for predictors; simple nonlinearity model
%   HRFtime     number of s to estimate HRF shape for, in FIR model.
%   delta       You must enter this if you want HRF shape estimate eff and you enter a model matrix X
%               Delta is n x m, with columns for conditions, made of ones and zeros, in TR-length time bins
%
% Required fields for condition function (vector) input:
%   ISI         sampling resolution of design vector, as time between elements in s
%   TR          sampling resolution of final design matrix (volumes acquired), in s
%
% Functions Called: (from OptimizeDesign GA toolbox)
% designvector2model.m
% tor_make_deconv_mtx2.m
% calcEfficiency.m

if ~isfield(P,'contrasts'), P.contrasts = eye(size(X,2));, end
if ~isfield(P,'HRFtime'), P.HRFtime = 30;, end
if isfield(P,'delta'), delta = P.delta;, end
if isfield(P,'TR'), TR = P.TR;, end
if isfield(P,'ISI'), ISI = P.ISI;, end

% -------------------------------------------------------------------------------------------------
% * Set up design matrix
% -------------------------------------------------------------------------------------------------

% For condition function
if any(size(X) == 1)    
    
    if ~isfield(P,'ISI'), error('ISI is required field in PARAMS input'), end
    if ~isfield(P,'TR'), error('TR is required field in PARAMS input'), end
    if ~isfield(P,'S'), P.S = [];, end
    if ~isfield(P,'HRF'), P.HRF = spm_hrf(.1);, HRF = HRF ./ max(HRF);, end
    if ~isfield(P,'nonlint'), P.nonlint = [];, end

    numsamps = ceil(length(X)*P.ISI/P.TR);
    
    if nargout > 1
        delta = [];
        for i = 1:max(double(X))
            delta(:,i) = (X == i);
        end
    end
    
    X = designvector2model(X,P.ISI,P.HRF,P.TR,numsamps,P.nonlint,P.S);
    
    if ~isfield(P,'Vi'), P.Vi = eye(size(X,1));, end
    

% For a design matrix
else                    
    if ~isfield(P,'S'), P.S = [];, end
    if ~isfield(P,'Vi'), P.Vi = eye(size(X,1));, end
    
end






if ~isempty(P.S), svi = P.S * P.Vi;, else, svi = P.Vi;, end

if ~isempty(P.contrasts)
        % -------------------------------------------------------------------------------------------------
		% * efficiency
		% -------------------------------------------------------------------------------------------------

        	xtxitx = pinv(X);                                       		% inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
            contrastweights = ones(1,size(P.contrasts,1));
			[dummy,varargout{1}] = calcEfficiency(contrastweights,P.contrasts,xtxitx,svi);
		   
else
        varargout{1} = [];
end
    
if isfield(P,'ISI') & isfield(P,'TR') & nargout > 1 & exist('delta') == 1
    
        % -------------------------------------------------------------------------------------------------
		% * HRF shape estimation efficiency
		% -------------------------------------------------------------------------------------------------

			[X2] = tor_make_deconv_mtx3(delta,round(P.HRFtime / P.TR),P.TR / P.ISI);
            if ~isempty(P.S), X2 = P.S * X2;,end
            
        	xtxitx = pinv(X2);                                       		% inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
			[dummy,varargout{2}] = calcEfficiency([],[],xtxitx,svi);

elseif nargout > 1
            disp(['Missing ISI or TR field, or already-constructed model entered: HRF efficiency not calculated.'])
            varargout{2} = [];
end

if nargout > 2
    varargout{3} = X; varargout{4} = X2;
    P.hrfsamples = round(P.HRFtime / P.TR);
    P.hrftimeres = P.TR / P.ISI;
    varargout{5} = P;
end

return
