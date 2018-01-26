function r = atlas2region(dat)
% Convert an atlas object to a region object
% r = atlas2region(dat)
%
% to-do: could add 'separate_into_contiguous' option
% would re-parse regions and adjust labels accordingly.

dat.dat = single(dat.dat); % int32 doesn't work, don't know why.

r = region(dat, 'unique_mask_values');


% Add labels to region names

if isempty(dat.labels)
    return
    
else
    
    for i = 1:length(r)
        
        if i > length(dat.labels)
            error('More regions than labels! Labels or integer codes may be wrong.');
        end
        
        r(i).shorttitle = dat.labels{i};
        r(1).custom_info1_descrip = dat.space_description;
        
        if ~isempty(dat.label_descriptions) && length(dat.label_descriptions) >= i
            
            r(i).title = sprintf('%s from %s', dat.label_descriptions{i}, dat.atlas_name);
            
        end
        
    end
    
end

end

