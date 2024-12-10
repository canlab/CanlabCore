function [param] = get_parameters_logit(hdrf,t,VM)
%
% [param] = get_parameters_logit(hdrf,t,VM)
%
% Estimate Height, time-to-peak and Width for the inverse logit model
%
% INPUT: hdrf,t,VM
% hdrf - HRF estimated using IL model
% t - vector of time points
% VM - IL model parameters 
%
% OUTPUT: param =[h,p,w]
% Height - h
% Time to peak - p (in time units of TR seconds)
% Width (at half peak) - w  
%
% By Martin Lindquist and Tor Wager
% Edited 12/12/06
%

% Estimate Heights and Time to peak:

[dL, dL2, d1, d2] = dilogit_dt2(t,VM);  % Calculate first and second derivativeof IL model
g = (abs(diff(sign(dL))) >= 1);
cnt = max(1,ceil(VM(2)));               %Peak has to be further along then T1 

p =-1; 
while(cnt < max(t)),
    if((g(cnt) == 1) & (dL2(cnt) <= 0)),
        p = cnt;
        h = hdrf(p);
        cnt = max(t) +1;
    end;
    cnt = cnt+1;
end;

if (p == -1),
    [h,p] = max(hdrf);
end;


% Interpolate to get finer estimation of p
if (p>1 & p<length(hdrf)),
    tp = p - 2 + (1:0.01:3);
    [h,delta] = max(Get_Logit(VM,tp));
    p = p-1 + delta/100;   
end;


% Estimate Width:

v = (hdrf >= h/2);
[a,b2] = min(diff(v));
v(b2+1:end) = 0;
[a,b1] = max(v);

% Interpolate to get finer estimation of v

if (p>1 && p<length(hdrf))
    tp = b1 - 1 + (0:0.01:1);
    htmp = Get_Logit(VM,tp);
    delta1 = sum((htmp >= h/2));
    tp = b2 + (0:0.01:1);
    htmp = Get_Logit(VM,tp);
    delta2 = sum((htmp >= h/2));
    w = sum(v) + delta1/100 + delta2/100;
else
    w = sum(v)+1;   
end;

% Set output vector param
param = zeros(3,1);
param(1) = h;
param(2) = p;
param(3) = w;

return;
