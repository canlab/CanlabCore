function delete_ana_imgs(imgs)
% Function to delete a list of .img files 
%   - automatically removes associated .hdr files as well.
%   - accepts cellstr or char inputs
%
% :Usage:
% ::
%
%     delete_ana_imgs(imgs)

    imgs = cellstr(imgs);
    for i = 1:length(imgs)
        delete(imgs{i});
        delete([imgs{i}(1:end-4) '.hdr']);
        if(exist([imgs{i}(1:end-4) '.mat'], 'file'))
            delete([imgs{i}(1:end-4) '.mat']);
        end
    end
end
