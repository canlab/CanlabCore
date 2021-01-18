
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [jd,pane]=msgError(txt, pause, where, title)
if nargin<4
    title='Error...';
    if nargin<3
        where='center';
        if nargin<2
            pause=0;
        end
    end
end
[jd, pane]=msg(txt, pause, where, title, 'error.png');
end
