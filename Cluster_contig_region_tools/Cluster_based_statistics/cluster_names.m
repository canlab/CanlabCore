function [cl,names] = cluster_names(cl,addflag, varargin)
% function [cl,names] = cluster_names(cl,[addflag], [overlay],[rename])
%
% Assign names to cl(x).shorttitle
% Do not use spaces or underscores or special chars for best results
%
% Addflag is optional; if 1, uses current orthview display
%
% If you have shorttitle field but want to rename, use rename flag = 1

spm_orthviews('Xhairs','on');

rename = [];
if nargin > 3, rename = varargin{2}; end

ovl = [];
if nargin > 2, ovl = varargin{1}; end

names = {};

dosurf = 0;

if nargin > 1 && addflag
    % don't create new figure
else
     cluster_orthviews(cl,{[1 0 0]}, 'solid', 'overlay', ovl);
end

for i = 1:length(cl)
    
    goname = 1;
    
    if dosurf, create_figure('Surface'); cluster_surf(cl(i), 5); end
    
    if isfield(cl(i),'shorttitle') && ~isempty(cl(i).shorttitle)
            goname = 0;
    end
    
    if rename
        goname = 1;
    end
        
    if goname
        spm_orthviews('Reposition',cl(i).mm_center);
        
        cl(i).shorttitle = input('Enter short name for this cluster: ','s');
        
        names{i} = cl(i).shorttitle;
    end
end


return