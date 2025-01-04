function [VM, hrf, fit, e, param] = Det_Logit(V0,t,tc,Run)
%
% [VM, h, fit, e, param] = Det_Logit_allstim(V0,t,tc,Run)
%
% Estimate inverse logit (IL) HRF model 
% Creates fitted curve - 3 logistic functions to be summed together - from parameter estimates
%
% INPUT: V0, t, tc, Run
% Run = stick function
% tc = time course
% t = vector of time points
% V0 = initial value for the parameter vector
%
% By Martin Lindquist, Christian Waugh and Tor Wager
% Created by Martin Lindquist on 10/02/09
% Last edited: 05/26/10 (ML)


numstim = length(Run);
len = length(Run{1});

% LB = [0.05, 1, 0, 0.05, 4, 0, 10];      % Lower bounds for parameters
% UB = [10, 15, 10, 10, 15, 5, 50];           % Upper bounds for parameters
% LB = repmat(LB, 1, numstim);
% UB = repmat(UB, 1, numstim);

% Remove intercept

b0 = pinv(ones(length(tc),1))*tc;
tc = tc - b0;

% Find optimal values

options = optimset('MaxFunEvals',10000000,'Maxiter',10000000,'TolX',1e-8,'TolFun',1e-8,'Display','off');

%VM = fminsearchbnd(@cost_allstim, V0, LB,UB,options,t,tc,Run);
VM = fminsearch(@msq_logit,V0,options,Run,t,tc);

% Use optimal values to fit hemodynamic response functions
hrf =zeros(length(t),numstim);
fitt = zeros(len,numstim);
param = zeros(3,numstim);

for g = 1:numstim
    hrf(:,g) = il_hdmf_tw2(t,VM(((g-1)*7+1):(g*7)));                   % Calculate HRF estimate (fit, given theta)
    param(:,g) = get_parameters2(hrf(:,g),t(end));
    fits(:,g) = conv(Run{g}, hrf(:,g));
    fitt(:,g) = fits(1:len,g);
end

fit = sum(fitt,2);
e = tc-fit;
fit = fit + b0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SUBFUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m=msq_logit(V,Run, t, tc)

numstim = length(Run);
len = length(Run{1});
h = zeros(length(t),numstim);
yhatt =zeros(len,numstim);

for k = 1:numstim
    h(:,k) = il_hdmf_tw2(t,V(((k-1)*7+1):(k*7)));           % Get IL model corresponding to parameters V
    yhat(:,k) = conv(Run{k}, h(:,k));                     % Convolve IL model with stick function
    yhatt(:,k) = yhat(1:len,k);
end

yhat2 = sum(yhatt,2); %Sum models together to get overall estimate

m = sum((tc-yhat2).^2);              % Calculate cost function

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h,base] = il_hdmf_tw2(t,V)
% inverse logit -- creates fitted curve from parameter estimates
%
% t = vector of time points
% V = parameters

% 3 logistic functions to be summed together
base = zeros(length(t),3);
A1 = V(1);
T1 = V(2);
d1 = V(3);
A2 = V(4);
T2 = V(5);
A3 = V(6);
T3 = V(7);
d2 = -d1*(ilogit(A1*(1-T1)) - ilogit(A3*(1-T3)))/(ilogit(A2*(1-T2)) + ilogit(A3*(1-T3)));
d3 = abs(d2)-abs(d1);

base(:,1)= d1*ilogit(A1*(t-T1))';
base(:,2)= d2*ilogit(A2*(t-T2))';
base(:,3)= d3*ilogit(A3*(t-T3))';
h = sum(base,2)';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L] = ilogit(t)
L = exp(t)./(1+exp(t));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
