function isok = check_extracted_data(cl)
% check_extracted_data Verify region-average data via re-extraction from source images.
%
% Re-extract the region-average data from 5 randomly chosen regions
% using spm_get_data.m and compare the result to the already-saved
% values in cl.dat. Useful as a sanity check when space/programming
% issues are suspected.
%
% :Usage:
% ::
%
%     isok = check_extracted_data(cl)
%
% :Inputs:
%
%   **cl:**
%        A valid region-class object (see region.m). cl(1).source_images
%        must still be on the path so the original data can be re-read.
%
% :Outputs:
%
%   **isok:**
%        Logical scalar: true if the correlation between re-extracted
%        and stored region averages is greater than 0.999 for every
%        sampled region.
%
% :Notes:
%
% You should not need to run this regularly, but you should if you
% suspect things have gone awry.
%
% :See also:
%   - region
%   - spm_get_data
%   - extract_roi_averages

isok = 1;

% Check that source images exist
imgs = check_valid_imagename(cl(1).source_images, 1);

wh_cl = randperm(length(cl));

for cl_num = wh_cl(1:5)
  
% Extract data from orig voxels
% ---------------------------------
Y = spm_get_data(imgs, cl(cl_num).XYZ);

avgdat = [nanmean(Y')' cl(cl_num).dat];

rr = corrcoef(avgdat);

isok = rr(1, 2) > .999;

okstr = {'No' 'Yes'};
fprintf('Checked data in region %3.0f, ok = %s\n', cl_num, okstr{1 + (rr(1, 2) > .999)});

end

end % function
