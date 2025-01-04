function [param] = get_parameters2(hdrf,t)
%
% Find model parameters
%
% Height - h
% Time to peak - p (in time units of TR seconds)
% Width (at half peak) - w  


% Calculate Heights and Time to peak:

% delta = 1/(t(2)-t(1));
% n = round(t(end)*0.6*delta)
n = round(length(t)*0.8);

[~,p] = max(abs(hdrf(1:n)));
h = hdrf(p);

%if (p > t(end)*0.6*delta), warning('Late time to peak'), end;
if (p > t(end)*0.8), warning('Late time to peak'), end;

if (h >0)
    v = (hdrf >= h/2);    
else
    v = (hdrf <= h/2);
end;
    
[~,b] = min(diff(v));
v(b+1:end) = 0;
w = sum(v);

cnt = p-1;
g =hdrf(2:end) - hdrf(1:(end-1));
while((cnt > 0) && (abs(g(cnt)) <0.001)),
    h = hdrf(cnt);
    p = cnt;
    cnt = cnt-1;
end;


param = zeros(3,1);
param(1) = h;
param(2) = p;
param(3) = w;

return;
