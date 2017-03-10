function [L] = ilogit(t)
% Calculate the inverse logit function corresponding to the value t
%
% :Usage:
% ::
%
%     function [L] = ilogit(t)
%
% :Output:
%
%   **L:**
%        exp(t)./(1+exp(t));
%
% ..
%    By Martin Lindquist and Tor Wager
%    Edited 12/12/06
% ..

L = exp(t)./(1+exp(t));
