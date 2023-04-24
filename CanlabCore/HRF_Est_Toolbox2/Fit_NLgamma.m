function [hrf, fit, e, param] =  Fit_NLgamma(tc, TR, Run, T )
%
% Fits non-linear Gamma HRF model  
%
% INPUTS:
% 
% tc    - time course
% TR    - time resolution
% Runs  - expermental design
% T     - length of estimated HRF
% 
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual time course
% param - estimated amplitude, height and width
%
% By Martin Lindquist
% Created by Martin Lindquist on 03/13/23


numstim = length(Run);
len = length(Run{1});
t=1:TR:T;

V0 = repmat([1 6 1]',1,numstim);  % [Height, Delay, Onset];

% Find optimal values

options = optimset('MaxFunEvals',10000,'Maxiter',10000,'TolX',1e-6,'TolFun',1e-6,'Display','off');
VM = fminsearch(@msq_nl_gamma,V0,options,Run,TR, T,tc);

% Use optimal values to fit hemodynamic response functions
hrf =zeros(length(t),numstim);
fitt = zeros(len,numstim);
param = zeros(3,numstim);

for g = 1:numstim
    hrf(:,g) = NL_gamma(TR, T, VM(:,g));                   % Calculate HRF estimate (fit, given theta)
    param(:,g) = get_parameters2(hrf(:,g),t);
    fits = conv(Run{g}, hrf(:,g));
    fitt(:,g) = fits(1:len);
end

fit = sum(fitt,2);
e = tc-fit;
%fit = fit + b0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SUBFUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m=msq_nl_gamma(V, Run, TR, T, tc)

numstim = length(Run);
len = length(Run{1});
t=1:TR:T;
h = zeros(length(t),numstim);
yhatt =zeros(len,numstim);

for k = 1:numstim
    h(:,k) = NL_gamma(TR, T, V(:,k));           % Get NL gamma model corresponding to parameters V
    yhat = conv(Run{k}, h(:,k));                     % Convolve IL model with stick function
    yhatt(:,k) = yhat(1:len);
end

yhat2 = sum(yhatt,2); %Sum models together to get overall estimate

m = sum((tc-yhat2).^2);              % Calculate cost function

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h] = NL_gamma(TR, T, V)
% inverse logit -- creates fitted curve from parameter estimates
%
% t = vector of time points
% V = parameters

% 3 logistic functions to be summed together

height = V(1);
%delay = 6;
delay = V(2);
udelay = 16;
%dispersion = 1;
dispersion = V(3);
udisp = 1;
rtou = 6;
onset = 1;
klength = T-1;

 %	p(1) - delay of response (relative to onset)	   6
 %	p(2) - delay of undershoot (relative to onset)    16
 %	p(3) - dispersion of response			   1
 %	p(4) - dispersion of undershoot			   1
 %	p(5) - ratio of response to undershoot		   6
 %	p(6) - onset (seconds)				   0
 %	p(7) - length of kernel (seconds)		  32
    
normhrf = spm_hrf(TR,[delay udelay dispersion udisp rtou onset klength]);
h = height.*normhrf ./ max(normhrf);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
