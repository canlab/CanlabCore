function [param, hdrf2] = Get_parameters_Logit2_2hrfs(hdrf,t,VM,t1t2)
%
% [param] = get_parameters_logit(hdrf,t,VM) 
%
% Different from version 1 in that it only calculates increases as 'peaks'
% Estimate Height, time-to-peak and Width for the inverse logit model
%
% INPUT: hdrf,t,VM (all with two rows for each HRF of interest)
% hdrf - HRF estimated using IL model
% t - vector of time points
% VM - IL model parameters 
% t1t2 - number of time points separating the onset of the first and second
% hrfs

% OUTPUT: param =[h,p,w]
% Height - h
% Time to peak - p (in time units of TR seconds)
% Width (at half peak) - w  
%
% By Martin Lindquist and Tor Wager
% Edited 12/12/06
%
% In addition, it'll calculate the width from two HRFs put together in case
% your design is such that two events are dependent on each other (i.e.
% cue->stimulus) 
%

% added by Christian Waugh

% Estimate Heights and Time to peak for first HRF:

[dL, dL2, d1, d2] = dilogit_dt2(t,VM(1, :));  % Calculate first and second derivativeof IL model
g = (abs(diff(sign(dL))) >= 1);
%cnt = max(1,ceil(VM(1, :)(2)));               %Peak has to be further along then T1 
cnt = max(1,ceil(min(min(VM(1,2),VM(1,5)),VM(1,7)))); 

p =-1; 

while(cnt < max(t)),
    if((g(cnt) == 1) && (dL2(cnt) <= 0)),
       p = cnt;
        h = hdrf(p,1);
        cnt = max(t) +1;
    end;
    cnt = cnt+1;
end;


if (p == -1) 
    [h,p] = max(hdrf(:,1));
end;


% Interpolate to get finer estimation of p
if (p>1 && p<length(hdrf(:,1))),
    tp = p - 2 + (1:0.01:3);
    [h,delta] = max(Get_Logit(VM(1, :),tp));
    p = p-1 + delta/100;   
end


% Estimate Width:


v = (hdrf(:,1) >= h/2);
[a,b2] = min(diff([v;0]));
if b2 < length(v)
v(b2+1:end) = 0;
end
[a,b1] = max(v);

% Interpolate to get finer estimation of v

if (p>1 && p<length(hdrf(:,1))),
    tp = b1 - 1 + (0:0.01:1);
    htmp = Get_Logit(VM(1, :),tp);
    delta1 = sum((htmp >= h/2));
    tp = b2 + (0:0.01:1);
    
    htmp = Get_Logit(VM(1, :),tp);
    delta2 = sum((htmp >= h/2));
    w = sum(v) + delta1/100 + delta2/100;
else
    w = sum(v)+1;   
end;

% Set output vector param
param = zeros(7,1);
param(1) = h;
param(2) = p;
param(3) = w;


% Estimate Heights and Time to peak for second HRF:

[dLb, dL2b, d1b, d2b] = dilogit_dt2(t,VM(2, :));  % Calculate first and second derivativeof IL model
gb = (abs(diff(sign(dLb))) >= 1);
%cnt = max(1,ceil(VM(1, :)(2)));               %Peak has to be further along then T1 
cntb = max(1,ceil(min(min(VM(2,2),VM(2,5)),VM(2,7)))); 

pb =-1; 
hb = 0;


while(cntb < max(t)),
    if((gb(cntb) == 1) && (dL2b(cntb) <= 0)),
        pb = cntb;
        hb = hdrf(pb,2);
        cntb = max(t) +1;
    end;
    cntb = cntb+1;
end;


if (pb == -1)
    [hb,pb] = max(hdrf(:,2));
end;


% Interpolate to get finer estimation of p
if (pb>1 & pb<length(hdrf(:,2))),
    tpb = pb - 2 + (1:0.01:3);
    [hb,delta] = max(Get_Logit(VM(2, :),tpb));
    pb = pb-1 + delta/100;   
end;


% Estimate Width:

vb = (hdrf(:,2) >= hb/2);

[ab,b2b] = min(diff([vb;0]));
if b2b < length(vb)
vb(b2b+1:end) = 0;
end
[ab,b1b] = max(vb);

% Interpolate to get finer estimation of v

if (pb>1 & pb<length(hdrf(:,2))),
    tpb = b1b - 1 + (0:0.01:1);
    htmpb = Get_Logit(VM(2, :),tpb);
    delta1 = sum((htmpb >= hb/2));
    tpb = b2b + (0:0.01:1);
    
    htmpb = Get_Logit(VM(2, :),tpb);
    delta2 = sum((htmpb >= hb/2));
    wb = sum(vb) + delta1/100 + delta2/100;
else,
    wb = sum(vb)+1;   
end;

% Set output vector param
param(4) = hb;
param(5) = pb;
param(6) = wb;


% Combine hrfs 1 and 2 to get an overall width estimate

hdrf2 = [hdrf(:,1); repmat(0,t1t2,1)] + [repmat(0,t1t2,1); hdrf(:,2)];

    vc = (hdrf2 >= h/2);

[ac,b2c] = min(diff([vc;0]));
if b2c < length(vc)
vc(b2c+1:end) = 0;
end
[ac,b1c] = max(vc);


wc = sum(vc)+1;   

param(7) = wc;

return;
