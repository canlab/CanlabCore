

% ALIGNK  calculate kernel alignment
%
% [X]=ALIGNK(K1,K2) calculates alignment between kernels K1 and K2.

function [x]=alignk(ik,ok)

x=sum(sum( ik .* ok)) / sqrt(sum(sum(ik.^2)) * sum(sum(ok.^2)) );