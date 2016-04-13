function varargout = efficency(X,P)
% function [con_eff,hrf_eff] = efficency(X,PARAMS)
% Tor Wager, 9/30/02
%
% a useful summary function for a design, incorporating tests of contrasts and of HRF shape estimation
%
% PARAMS structure contains inputs
% Optional fields:
%   contrasts   required for testing contrast efficiency; see calcEfficiency.m
%   Vi          intrinsic autocorrelation matrix; empty assumes identity matrix (independence)
%   S           smoothing matrix; empty for no smoothing
%   HRF         canonical hemodynamic response function; empty assumes SPM's canonical HRF
%               in this function, FIR estimation length fixed at 12 s.
%   nonlint     saturation threshold (ceiling) for predictors; simple nonlinearity model
%
% Required fields for condition function (vector) input:
%   ISI         sampling resolution of design vector, as time between elements in s
%   TR          sampling resolution of final design matrix (volumes acquired), in s
%
% Functions Called: (from OptimizeDesign GA toolbox)
% designvector2model.m
% tor_make_deconv_mtx2.m
% calcEfficiency.m

hrfsec = 30;

if ~isfield(P,'contrasts'), P.contrasts = []; end

if ~isfield(P,'S'), P.S = []; end

if any(size(X) == 1)    % then it's a condition function
    
    if ~isfield(P,'ISI'), error('ISI is required field in PARAMS input'), end
    if ~isfield(P,'TR'), error('TR is required field in PARAMS input'), end
    if ~isfield(P,'Vi'), P.Vi = eye(size(des,1)); end
    
    if ~isfield(P,'HRF'), P.HRF = spm_hrf(.1); end
    if ~isfield(P,'nonlint'), P.nonlint = []; end

    numsamps = ceil(length(X)*P.ISI/P.TR);
    X = designvector2model(X,P.ISI,P.HRF,P.TR,numsamps,P.nonlint,P.S);
    
    delta = [];
    for i = 1:max(X)
        delta(:,i) = (X == i);
    end
else                    % then it's a design matrix
end
    


svi = P.S * P.Vi;

if ~isempty(P.contrasts)
        % -------------------------------------------------------------------------------------------------
		% * efficiency
		% -------------------------------------------------------------------------------------------------

        	xtxitx = pinv(X);                                       		% inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
            
            contrastweights = ones(1,size(xtxitx,1)-1);
			[dummy,varargout{1}] = calcEfficiency(contrastweights,P.contrasts,pinv(X),svi);
		   
else

        [eff, eff_vector] = calcEfficiency([], [], pinv(X), svi);
        varargout{1} = eff_vector;
end
    
if isfield(P,'ISI') & isfield(P,'TR') & nargout > 1 & exist('delta', 'var')
    
        % -------------------------------------------------------------------------------------------------
		% * HRF shape estimation efficiency
		% -------------------------------------------------------------------------------------------------
        

      
			[X2] = tor_make_deconv_mtx2(delta,round(hrfsec / P.TR),P.TR / P.ISI);
            if ~isempty(P.S), X2 = P.S * X2;,end
            
        	xtxitx = pinv(X2);                                       		% inv(X'S'SX)*(SX)'; pseudoinv of (S*X)
			[dummy,varargout{2}] = calcEfficiency([],[],xtxitx,svi);

elseif nargout > 1
            disp(['Missing ISI or TR field, or already-constructed model entered: HRF efficiency not calculated.'])
            varargout{2} = [];
end

return