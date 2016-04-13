function [w,e] = noise_arp(n, phi, sigma, e)
% :Usage:
% ::
%
%     w = noise_arp(n, phi)
%
% Generate n-length AR(p) noise vector w
% given phi (default = [.5 .1])
%
% Based on Gaussian noise process e with standard deviation sigma (default = 1)
%
% w will have a variance greater than that of the underlyling gaussian
% process e
%
% ..
%    Martin Lindquist / Tor Wager, Feb 2007
% ..

    if nargin < 2
        phi = [0.5 .1];
    end

    if nargin < 3
        sigma = 1;
    end

    % number of AR parameters (model order):
    p = length(phi);

    if nargin < 4
        % underlying gaussian noise process:
        e = normrnd(0,sigma,n,1);
    end

    % w is colored noise:
    w = e;

    for i= (p+1):n
        arpart = phi * w(i-1 : -1 : i-p);

        w(i) = arpart + e(i);
    end

end
