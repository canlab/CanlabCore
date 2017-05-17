% function enlarge_axes(figure_handle, scale_factor, [zoom value])
%
% Shrink:
% enlarge_axes(gcf, .8)
%
% Grow:
% enlarge_axes(gcf, 1.3)
%
% Grow and zoom:
% enlarge_axes(gcf, 1.3, 1.2)
%
% Just zoom in:
% enlarge_axes(gcf, 1, 1.2)
%
function enlarge_axes(fh, scale_factor, varargin)
    
    zoomval = 0; 
    if length(varargin) > 0, zoomval = varargin{1}; end
    
    figure(fh);

    if nargin < 2, scale_factor = 1.35; end

    h = findobj(get(gcf, 'Children'), 'Type', 'axes');
    for i = 1:length(h)
        pp = get(h(i), 'Position');
%         
%         
%         set(h(i), 'OuterPosition', pp);
%         set(h(i), 'Position', pp);


        set(h(i), 'Position', [pp(1) pp(2) pp(3)*scale_factor pp(4)*scale_factor]);
        
        if zoomval
            axes(h(i))
            camzoom(zoomval)
        end
        
    end
end
