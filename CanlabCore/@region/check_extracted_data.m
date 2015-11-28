function isok = check_extracted_data(cl)
% Checks the data, just in case of space/programming issues, 
% by re-extracting the region average data from 5 random regions 
% using spm_get_data.m, and compares it to the already-saved values
%
%:Inputs:
%
%   **cl:**
%        must be a valid region object (see region.m)
%        and cl(1).source_images must still be on the path.
%
% You should not need to run this regularly -- but you should if you
% suspect things have gone awry.

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
