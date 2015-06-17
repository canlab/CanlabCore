function [X,d,out] = rt2delta(ons,rt,TR,varargin)
% [X,d,out] = rt2delta(ons,rt,TR,varargin)
%
% tor wager, 3/6/04
%
% given onset times in s and rt times in ms, and a TR rep. time for
% sampling, this function builds a model (X) that contains linear
% predictors convolved with a canonical hrf (SPM) for each trial type
% (trial types are stored as cells in the ons and rt cell arrays, one cell
% per trial type).  The model (X) also contains linear and quadratic predictors 
% for the effects of log(rt) on activation in each trial type.  d is an 
% indicator matrix of onset times in TRs for each event of interest, scaled
% by RT in the case of the linear and quadratic predictors.  The order in X
% and d is Activation-type1 Activation-type2 ... Act-type n LinearRT-Type
% 1...Linear-type n Quadratic-type 1 ... Quad type n.
%
% A third output argument contains separate matrices for each X and d for
% activation, linear RT, and quadratic RT effects, as well as a d and X for
% classification of trials into low, medium, and high RTs for each trial
% type, in that order.  This is an alternative to fitting linear and
% quadratic trends.
%
% Problems and solutions:
% - RT is typically skewed.  This function uses the Z(log RT) transform,
% and windsorizes Z values to 3 standard deviations, to avoid outliers.
% Squared RTs are additionally windsorized in a separate step.
%
% - RT may produce latency differences as well as magnitude of activation
% differences.  No solution in this algorithm
%
% - several responses may occur in the same TR.  This function builds an
% accurate d up to 2 responses per TR, although the 2nd (?) RT only is
% used.  This is a limitation, but most designs will hopefully avoid
% multiple trials in the same TR.


% kludgy fix for when you enter onsets in TRs instead of s, but want
% convolution to be right -- add 2nd "conv TR"
convTR = TR;
if length(varargin) > 0, convTR = varargin{1};, end

hrf = spm_hrf(convTR) ./ max(spm_hrf(convTR));

[X,d] = onsets2delta(ons,TR,ceil(max(cat(1,ons{:}))));  % fix length to longest
dout = zeros(size(d));
rtclass = [];

if length(rt) ~= length(ons), error('RT and onsets must be same length'),end

for i = 1:length(rt)
    
    dd = 1+round(ons{i} ./ TR); %find(d{i});    % index of onsets
    %dd = 1+round(ons{i} ./ 1);                   % keep at 1 to preserve orig input
    
    if length(dd) ~= length(rt{i}), warning(['Number of RTs for condition ' num2str(i) ' does not match onsets']),keyboard;end
    
    % replace any RT less than 100 ms with mean RT
    rts = rt{i};
    wh = find(rts < 100); 
    rts(wh) = mean(rts);
    
    % center and windsorize log-transformed rts
    rts = scale(log(rts),1);                     % log transform and center
    rts2 = rts .^ 2;                            % save squared term for later
    wh = find(abs(scale(rts)) > 3);
    rts(wh) = sign(rts(wh)) .* 3;
    
    % save linear term for output
    dout(dd,i) = rts;
    
    % define quadratic term (centered)
    rts2 = scale(rts2,1);   % center
    wh = find(abs(scale(rts2)) > 3);
    rts2(wh) = sign(rts2(wh)) .* 3;
    dquad(dd,i) = rts2;
        
    if nargout > 2, 
        % classify RTs into bins, bottom middle and top 3rd of Z-scores
        rts = scale(rts);   % z-score
        rtclass(:,end+1:end+3) = zeros(size(X,1),3);
        tmp = [rts dd]; tmp = sortrows(tmp,1);  % sort by RT time
        b1 = round(size(tmp,1) ./ 3);
        rtclass(tmp(1:b1,2),end-2) = 1;
        rtclass(tmp(b1+1:2*b1,2),end-1) = 1;
        rtclass(tmp(2*b1+1:end,2),end) = 1;
        %rtclass(dd(rts < -.4125),end-2) = 1;
        %rtclass(dd(rts >= -.4125 & rts <= .4125),end-1) = 1;
        %rtclass(dd(rts < .4125),end) = 1;
    end
          
end

if nargout > 2,  
    out.basicX = X; 
    out.basicd = d;, 
    out.rtclass = rtclass;, 
    out.rtclassX = getPredictors(rtclass,hrf);
    out.rtlinearX = getPredictors(dout,hrf);
    out.rtlineard = dout;
    out.rtquadX = getPredictors(dquad,hrf);
    out.rtquadd = dquad;
end

d = [cell2mat(d) dout dquad];
X = getPredictors(d,hrf);
X(:,end+1) = 1;

if nargout > 2,  out.rtX = X; out.rtd = d;, end



return

