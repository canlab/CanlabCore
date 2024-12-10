function [VM, h, fit, e, param] = Det_Logit(V0,t,tc,Run)
%
% [theta,HH,C,P] = Anneal_Logit(theta0,t,tc,Run)
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
% By Martin Lindquist and Tor Wager
% Edited 10/01/09
%

% Find optimal values
options = optimset('MaxFunEvals',10000000,'Maxiter',10000000,'TolX',1e-8,'TolFun',1e-8,'Display','off');
VM = fminsearch(@msq_logit,V0,options,Run,t,tc);
VM

% Use optimal values to fit hemodynamic response functions
h = il_hdmf_tw2(t,VM(1:7));

%[param] = get_parameters2(h,t);
[param] = get_parameters_logit(h,t,VM(1:7));

len = length(Run);
fit = conv(Run, h);
fit = fit(1:len);

e = tc-fit;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SUBFUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m=msq_logit(V,Run, t, tc)

HR = il_hdmf_tw2(t,V(1:7));
len = length(Run);
timecourse = conv(Run, HR);
timecourse = timecourse(1:len);

m=sum((tc-timecourse).^2);

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
% 
% function [param] = get_parameters2(hdrf,t)
% % Find model parameters
% %
% % Height - h
% % Time to peak - p (in time units of TR seconds)
% % Width (at half peak) - w  
% 
% % Calculate Heights and Time to peak:
% 
% n = ceil(t(end)*0.8);
% [h,p] = max(abs(hdrf(1:n)));
% h = hdrf(p);
% 
% if (h >0)
%     v = (hdrf >= h/2);    
% else
%     v = (hdrf <= h/2);
% end;
%     
% [a,b] = min(diff(v));
% v(b+1:end) = 0;
% w = sum(v);
% 
% cnt = p-1;
% g =hdrf(2:end) - hdrf(1:(end-1));
% while((cnt > 0) & (abs(g(cnt)) <0.001)),
%     h = hdrf(cnt);
%     p = cnt;
%     cnt = cnt-1;
% end;
% 
% param = zeros(3,1);
% param(1) = h;
% param(2) = p;
% param(3) = w;
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
