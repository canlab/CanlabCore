function [dLdt, d2Ldt2, d1, d2] = dilogit_dt2(t,V)
%
% [dLdt, d2Ldt2, d1, d2] = dilogit_dt2(t,V)
%
% Calculate first and second temoral derivative of inverse logit (IL) HRF model 
% curve
%
% INPUT: V, t
% t = vector of time points
% V = parameters
%
% OUTPUT: dLdt, d2Ldt2
% dLdt = first temporal derivative of IL model
% d2Ldt2 = second temporal derivative of IL model
%
% By Martin Lindquist and Tor Wager
% Edited 12/12/06
%

% Set parameter values
A1 = V(1);
T1 = V(2);
d1 = V(3);
A2 = V(4);
T2 = V(5);
A3 = V(6);
T3 = V(7);
d2 = -d1*(ilogit(A1*(1-T1)) - ilogit(A3*(1-T3)))/(ilogit(A2*(1-T2)) + ilogit(A3*(1-T3)));
d3 = abs(d2)-abs(d1);

% Calculate the first and second derivative
dLdt = d1*A1*ilogit(A1*(t-T1))./(1+exp(A1*(t-T1))) + d2*A2*ilogit(A2*(t-T2))./(1+exp(A2*(t-T2))) + d3*A3*ilogit(A3*(t-T3))./(1+exp(A3*(t-T3)));
d2Ldt2 = d1*(A1^2)*(exp(A1*(t-T1)) - exp(2*A1*(t-T1)))./((1+exp(A1*(t-T1))).^3) + d2*(A2^2)*(exp(A2*(t-T2)) - exp(2*A2*(t-T2)))./((1+exp(A2*(t-T2))).^3) + d3*(A3^2)*(exp(A3*(t-T3)) - exp(2*A3*(t-T3)))./((1+exp(A3*(t-T3))).^3); 

return;