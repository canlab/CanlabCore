function check_roi_extraction()

% test ROI extraction
mask_image = fmri_data(which('atlas_labels_combined.img'));
dat = mask_image;
wh_region = dat.dat;

% create timeseries data -- add data "after" 1st image
nimgs = 50;

n = unique(wh_region);
ts = linspace(-10, 10, nimgs);
v = size(dat.dat, 1);
dat.dat = zeros(v, nimgs);

for i = 1:length(n)
    my_indx = wh_region == n(i);

    dat.dat(my_indx, :) = repmat(ts .* sqrt(i), sum(my_indx), 1);
end

fprintf('\nDATA GENERATED\n\n');

cl1 = extract_roi_averages(dat, mask_image, 'unique_mask_values');

% extracted average values for all the regions in mask_image. 
all_reg1 = cat(2, cl1(:).dat);

dat.fullpath = fullfile(pwd, 'test_image.img');
write(dat);
fprintf('\nSAVED DATA TO %s\n\n', dat.fullpath);

% reload
fprintf('\nLOADING DATA FROM %s\n\n', dat.fullpath);
dat = fmri_data(dat.fullpath, which('gray_matter_mask.img'));

% region obj -> extract_data
% test extract_roi_averages
cl2 = extract_roi_averages(dat, mask_image, 'unique_mask_values');

% extracted average values for all the regions in mask_image. 
all_reg2 = cat(2, cl2(:).dat);

% uniq_reg = sort(unique(mask_image.dat));

% compare roi averages of initially data generated and data loaded from file
if is_equal(all_reg1, all_reg2, 0.001)
    fprintf('\nExtracted ROI averages are equal with an error threshold of 0.001\nPASS\n');
else
    fprintf('\nExtracted ROI averages are not equal with an error threshold of 0.001\nFAIL\n');
end


% ------------------------------------------------------
% INLINE FUNCTION
% ------------------------------------------------------

    function res = is_equal(A, B, error_threshold)
        err = abs((A-B)./B) <= error_threshold;

        res = all(err(:));
    end % is_equal

end  % main function