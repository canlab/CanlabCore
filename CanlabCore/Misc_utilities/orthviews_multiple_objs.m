function orthviews_multiple_objs(imgs)
% plot multiple image objects on one orthviews
%
% :Input:
%
%   cell array of image_vectors or statistic_images
%
% ..
%    Yoni Ashar, 11/2014
% ..

    n=length(imgs);

    overlay = which('SPM8_colin27T1_seg.img');

    spm_check_registration(repmat(overlay, n, 1));

    handle_indices = 1:n;

    for i=1:n
        orthviews(imgs{i}, 'han',i);
        %out{i} = o{i};
    end

end
