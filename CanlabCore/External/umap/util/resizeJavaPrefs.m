%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [d, resized]=resizeJavaPrefs(app, J, H)
d=J.getPreferredSize;
if app.toolBarFactor>0
    d=java.awt.Dimension(d.width*.55, d.height*.6);
    J.setPreferredSize(d);
    resized=true;
else
    resized=false;
end
end
