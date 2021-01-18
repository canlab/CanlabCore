%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function pu=reportClusterMerging(loop, changes, peaks, pu, startTime)
cur=get(0, 'CurrentFigure');
if loop==1 && peaks>200
    pu=PopUp(...
        sprintf('<html>Merging valleys between %d cluster peaks', peaks), ...
        'south east', 'Cluster valley merging');
    if ~isempty(cur)
        figure(cur);
    end
    return;
end
m=mod(loop,100);
if m==1
    str='st';
elseif m==2
    str='nd';
elseif m==3
    str='rd';
else
    str='th';
end

txt=sprintf('%s changes require %d%s merging loop\n', ...
    String.encodeInteger(changes), loop, str);
if isempty(pu)
    timeSpentSoFarClustering=toc(startTime);
    if timeSpentSoFarClustering>2 % more than 2 seconds??
        if changes>1000 || (loop>2&&changes>400) || loop>4
            pu=PopUp(txt,'south east', 'Cluster valley merging');
            if ~isempty(cur)
                figure(cur);
            end
        end
    end
elseif ~isempty(pu)
    pu.setText(txt);
    if ~isempty(cur)
        figure(cur);
    end
end
if Density.IsDebugging
    fprintf(txt);
end
end