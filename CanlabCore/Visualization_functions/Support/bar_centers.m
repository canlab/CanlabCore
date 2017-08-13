% function centers = bar_centers(handles)
%
% BAR_CENTERS determines the actual X centers of grouped bars for further plotting purposes
%   Inputs:
%       handles - a handle or handles to start searching from. Ideally, this should be the barseries handle(s), but axes handles will usually work, too
%   Outputs:
%       centers - a list of centers, one row per group

function centers = bar_centers(handles)

    hPatches = findobj(handles, 'Type', 'patch'); % old Matlab, pre-2015
    
    if isempty(hPatches)
       % new matlab?
       hPatches = findobj(handles, 'Type', 'bar');
    end
    
    for i=1:length(hPatches)
        centers(i,:) = mean(get(hPatches(i), 'XData'));
    end
    centers = flipud(centers);
    
end