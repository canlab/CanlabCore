%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Funded by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function val = clip(val)
%CLIP Standard clamping of a value into a fixed range (in this case -4 to
% 4). For speed reasons, this function is no longer used by the MATLAB
% algorithm.
%
% val = CLIP(val)
% 
% Parameters
% ----------
% val: double
%     The value to be clamped.
% 
% Returns
% -------
% The clamped value, now fixed to be in the range -4 to 4.

 val = min(4, val);
 val = max(-4, val);