function orthviews_multiple_objs(imgs, varargin)
% plot multiple image objects on one orthviews
%
% :Input:
%
%   cell array of image_vectors or statistic_images. varargin gets passed
%   to orthviews -- see options there
%
% ..
%    Yoni Ashar, 11/2014
% ..

    n=length(imgs);

    overlay = which('fmriprep20_template.nii.gz');

    spm_check_registration(repmat(overlay, n, 1));

    handle_indices = 1:n;

    for i=1:n
        orthviews(imgs{i}, 'han',i, varargin{:});
        %out{i} = o{i};
    end

end
