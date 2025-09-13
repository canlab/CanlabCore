function [V, K, acf_true] = arp_V(rho, n)
% Generate autocorrelation matrix V, whitening matrix K, autocorrelation
% function acf_true from a series of AR(p) coefficients rho
%
% rho = vector of AR(p) coefficients
% n = length of time series (num observations)
%
% Example:
% [V, K, acf_true] = arp_V([.7 .1], 100)
% NOTE: IN DEVELOPMENT, DO NOT USE WITHOUT FURTHER CHECKING/WORK
%
% [V, K, acf_true] = arp_V([.6 .2], 100); create_figure('V'); plot(acf_true); set(gca, 'XLim', [1 10]);
%
% SPM version:
% [Q] = spm_Q([.6 .2], 100); 
% acf_true = Q(1, :); 
% create_figure('V'); plot(acf_true); set(gca, 'XLim', [1 10]);

% Programmers' Notes:
% by Tor Wager, with input from Martin Lindquist
% -----------------------------------------------
% Create autocorrelation matrix V
% This is same as spm_Q., but i'm not sure spm_Q works right for ar_p ? 
% e(i+1) = rho * e(i) + eta(i), e = error, eta = innovation/new random error

acf_true = zeros(1, n);
acf_true(1) = 1;

for rho_p = 1:length(rho) % for each ar coefficient

    for i = (rho_p + 1):n

        acf_true(i) = acf_true(i) + rho(rho_p) .^ abs(i - rho_p);

    end

end % ar coeff

V = toeplitz(acf_true);

K = inv(V .^ 0.5);
% -----------------------------------------------

end