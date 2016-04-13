function dice_coeff = dice_coeff_image(imagelist)
% dice coefficient for the overlap between each pair of a set of images
%
% :Usage:
% ::
%
%     dice_coeff = dice_coeff_image([imagelist or fmri_data object])
%
% ..
% tor wager, June 2010
% ..

if ischar(imagelist) % Read data
    obj = fmri_data(imagelist);
    imagelist = obj;
    
elseif isa(imagelist, 'fmri_data') || isa(imagelist, 'statistic_image') || isa(imagelist, 'image_vector')
    
    % do nothing
    
else
    error('Unknown data type for imagelist.  Enter char array of filenames or image object.');
    
end

N = size(imagelist.dat, 2);
if N == 1, error('Must enter at least two images'); end

%     [maskInfo, dat] = iimg_read_img(imagelist(1, :));
%
%     for i = 2:N
%
%         dat2 = scn_map_image(imagelist(i, :), imagelist(1, :));
%         dat2 = dat2(:);
%
%         dat = [dat dat2];
%
%     end


dat = imagelist.dat;

dat = double(abs(dat) > eps);

n = dat'*dat; % intersection * 2

% simple, but not very efficient
for i = 1:N
    for j = 1:N
        
        d(i, j) = sum(dat(:, i) | dat(:, j)); % union
        
    end
end

dice_coeff = n ./ d;

end

