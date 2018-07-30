function zoom_in_on_regions(o2, cl, orientation)
% Assumes you have one axis per region in a montage registered in an
% fmridisplay object (o2). Given the object o2, which montage to adjust
% (montage_num), a corresponding region object with coordinates (cl), and an
% orientation ('axial', 'saggital', or 'coronal'), this object method zooms 
% in on each cluster by adjusting the axes.  
%
% Tor Wager, Feb 2018
%
% Examples:
% thal = load_atlas('thalamus');
% o2 = montage(thal, 'nosymmetric', 'regioncenters');
% zoom_in_on_regions(o2, 1, atlas2region(thal), 'axial');

% Zoom in on regions!

border_prop = .9;  % proportion border

for i = 1:length(cl)
    
    xyzminmax = [min(cl(i).XYZmm, [], 2) max(cl(i).XYZmm, [], 2)];  % 3 x 2, 3 min  3 max mm coords
    
    switch orientation
        case 'axial'
            xl = xyzminmax(1, :);
            yl = xyzminmax(2, :);
            
        case 'saggital'
            xl = xyzminmax(2, :);
            yl = xyzminmax(3, :);
            
        case 'coronal'
            xl = xyzminmax(1, :);
            yl = xyzminmax(3, :);
            
    end
    
    myborder = border_prop * range(xl);     % multiple of range from min to max coordinate
    myborder = max([myborder 20]);          % border at least 10 mm
    xl(1) = xl(1) - myborder;
    xl(2) = xl(2) + myborder;
    
    myborder = border_prop * range(yl);
    myborder = max([myborder 20]);          % border at least 10 mm
    yl(1) = yl(1) - myborder;
    yl(2) = yl(2) + myborder;
    
    %my_handle = o2.montage{montage_num}.axis_handles(i);
    my_handle = o2.montage{i}.axis_handles;
    
    % fix in case lower and upper are the same
    % this should no longer happen with min border of x mm
%     if ~diff(xl), xl = [xl(1) - 5 xl(1) + 5]; end 
%     if ~diff(yl), yl = [yl(1) - 5 yl(1) + 5]; end 
    
    set(my_handle, 'XLim', xl, 'YLim', yl);     % which of the 3 vals should be chosen depends on slice orientation (sagg, cor, ax)
    
end

end % function

