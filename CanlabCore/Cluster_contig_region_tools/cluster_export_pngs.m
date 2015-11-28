function cl = cluster_export_pngs(cl,useexisting,overlay,xhairson)
% Save png images of SPM orthviews windows for each cluster in a set (cl
% structure)
%
% :Usage:
% ::
%
%    cl = cluster_export_pngs(cl,[useexisting],[overlayimagename],[xhairson])
%
% names from cl(x).shorttitle are used
% useexisting is optional: 1 uses existing orthviews display (default), 0 creates a
% new one with the clusters
%
% :Example:
% ::
%
%    % use existing
%    cluster_export_pngs(cl, 1, EXPT.overlay);
%
% ..
%    tor wager, aug 3, 06
% ..

if nargin < 2, useexisting = 1; end
if nargin < 3, overlay = []; end
if nargin < 4, xhairson = 0; end

if ~useexisting
    cluster_orthviews(cl,'bivalent');
end

% subdirectory: prompt
dosubdir = input('Create a subdirectory for png images? (type name or return to use current dir) ','s');
if ~isempty(dosubdir)
    mkdir(dosubdir)
    cd(dosubdir)
end

% Make sure we have names or ask to create them
if ~isfield(cl(1),'shorttitle') || isempty(cl(1).shorttitle)
    donames = input('Name these clusters before saving? (1/0) ');
    if donames
        cl = cluster_names(cl,1);
    else
        for i = 1:length(cl)
            cl(i).shorttitle = [];
        end
    end
end

% set up orthviews figure
scn_export_spm_window('setup',overlay);

if xhairson, spm_orthviews('Xhairs','on'); end

% save png images of each cluster

for i = 1:length(cl)
    spm_orthviews('Reposition',cl(i).mm_center);
    spm_orthviews_showposition;
    
    % evaluate any existing window button up function (e.g., update
    % scatterplots)
    fh = findobj('Type','Figure','Tag','Graphics');
    funh = get(fh,'WindowButtonUpFcn');
    funh()

    scn_export_spm_window('save',['Cl' num2str(i) '_' cl(i).shorttitle]);
end

if ~isempty(dosubdir)
    cd ..
end

return
