function riverplot_remove_ribbons(ribbons)
% riverplot_remove_ribbons(ribbons)

for i = 1:size(ribbons, 1)
    for j = 1:size(ribbons, 2)
        
        if ~isempty(ribbons{i, j})
            delete(ribbons{i, j}.line1.h);
            delete(ribbons{i, j}.line2.h);
            delete(ribbons{i, j}.patchh);
        end
    end
end

end % function 