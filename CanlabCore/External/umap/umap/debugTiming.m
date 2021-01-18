%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
function debugTiming(description)
DEBUGGING=false;
if DEBUGGING && ~isdeployed
    msg([description ' ' num2str(toc)], 40, ...
        'south west+', 'Timing test', 'none')
end
end