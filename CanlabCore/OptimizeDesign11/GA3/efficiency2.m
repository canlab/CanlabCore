function [e,e2,ediff,esum] = efficiency2(x,tr,hp,rho)
% [e,e2,ediff,esum] = efficiency2(x,tr,hp,rho)
%
% hp = hp filter
% rho = AR(1) coefficient
%
% special for block_filter_sim.m

% scaling factor: make efficiency scale invariant
% assume intercept is the last predictor

n = size(x,1);


c = [1 -1 0]';

e = []; e2 = 1;
% manual filtering way, no fancy formulae
% for i = 1:2, X2(:,i) = hpfilter(x(:,i),tr,hp,n);, end
% X2(:,end+1) = 1;
% xtxi = inv(X2'*X2);
% e = 1./(sqrt(diag(xtxi)));
% e2 = 1./(sqrt(diag(c' * xtxi * c)));

% adjusts for scaling way
%p = [std(x(:,1:end-1)) sum(x(:,1))./n ];
%e2 = 1./(p'.*sqrt(diag(inv(x'*x))));



[S,V,svi,KH] = getSmoothing(hp,0,tr,n,[1 rho rho^2 rho^4]);

%V=eye(size(x,1)); 
x = S * x;
xtxi = inv(x'*x);
px = xtxi * x';

% var/covar
vcv = px * S * V' * V * S * px';

%ediff = 1./sqrt(diag(c' * px * V' * V * px' * c));
ediff = 1./sqrt(diag(c' * vcv * c));

%ediff = ediff ./ std(c);  % abs(c(1)) ??

c = [1 1 0]';
esum = 1./sqrt(diag(c' * vcv * c));
%esum = esum ./ std(c);

% the thing is that it takes var(predictors) out of the equation.
% predictor variance may really mean something in terms of power --
% but you have to be able to equate the scales for diff. predictors in some
% way.