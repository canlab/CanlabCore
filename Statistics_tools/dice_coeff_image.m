function dice_coeff = dice_coeff_image(imagelist)
    % dice_coeff = dice_coeff_image(imagelist)
    %
    % dice coefficient for the overlap between each pair of a set of images
    % tor wager, June 2010

    N = size(imagelist, 1);
    if N == 1, error('Must enter at least two images'); end

    [maskInfo, dat] = iimg_read_img(imagelist(1, :));

    for i = 2:N

        dat2 = scn_map_image(imagelist(i, :), imagelist(1, :));
        dat2 = dat2(:);

        dat = [dat dat2];

    end

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

