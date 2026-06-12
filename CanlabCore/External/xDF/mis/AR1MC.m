function [V,Z,P,R2Zcrt,arone]=AR1MC(Y,T)
% This function is AR(1) Monte-Carlo estimation of unbiased variance. 
% Copy-pasted from FSLnets toolbox. 
%
% NB! This function was designed for estimation of *global* correction (subjects)
% factors.
%
%%%REFERENCES:
%  Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. 
% "Effective Degrees of Freedom of the Pearson's Correlation Coefficient 
%  under Serial Correlation." bioRxiv (2018): 453795.
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

if  size(Y,1) ~= T
    Y = Y';
    %warning('Oi!')
end

nn = size(Y,2);

Y  = Y - mean(Y); %demean
%Estimate AR1--------------------------------------------------------------

arone0=[];
for i=1:nn
  g     = Y(:,i);  
  arone0 = [arone0 sum(g(1:end-1).*g(2:end))/sum(g.*g)];
end
arone  = arone0;
arone0 = median(arone0);

netmat = corr(Y);
clear Y
% create null data using the estimated AR(1) coefficient-------------------
for i=1:50 
  Y(1)=randn(1);
  for t=2:T %SA: creat random numbers with ar1 element for rest of the *run*
    Y(t)=Y(t-1)*arone0+randn(1);
  end
  Yr(:,i)=Y;
end

IDX     = find(triu(ones(50),1));
Yrc     = corr(Yr);

V       = var(Yrc(IDX));
%This is 1/N-3; where N-3 is estimated by the std of the corrs 
R2Zcrt  = 1./std(atanh(Yrc(IDX)));

Z       = atanh(netmat).*R2Zcrt;
P       = 2 * normcdf(-abs(Z)); 

