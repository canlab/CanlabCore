%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function ok=initJava
ok=false;
try
    edu.stanford.facs.swing.StochasticGradientDescent.EPOCH_REPORTS;
catch
    jar=fullfile(fileparts(mfilename('fullpath')), 'umap.jar');
    javaaddpath(jar);
end
try
    edu.stanford.facs.swing.StochasticGradientDescent.EPOCH_REPORTS;
    ok=true;
catch
end