%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

function screens=javaScreens
screens={};
ge=java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment;
physicalScreens=ge.getScreenDevices;
N=length(physicalScreens);
for i=1:N
    pe=physicalScreens(i).getDefaultConfiguration.getBounds;
    screens{end+1}=pe;
end
end
