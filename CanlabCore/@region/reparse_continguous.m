function clout = reparse_continguous(cl)
% reparse_continguous Re-define regions in a region object based on contiguous blobs.
%
% Convert the input region object to an image_vector representation,
% then re-build a region-class object array using contiguous-region
% parsing. If the input regions carry .all_data, the per-voxel data is
% redistributed across the new contiguous regions and re-averaged into
% .dat for each new region.
%
% :Usage:
% ::
%
%     clout = reparse_continguous(cl)
%
% :Inputs:
%
%   **cl:**
%        A region-class object array, typically resulting from
%        operations that may have produced regions with non-contiguous
%        voxels (e.g., posneg_separate).
%
% :Outputs:
%
%   **clout:**
%        A new region-class object array with one element per
%        contiguous blob, with .all_data and .dat re-distributed where
%        possible.
%
% :See also:
%   - region2imagevec
%   - posneg_separate
%   - region
%
% ..
%    NEEDS SOME ADDITIONAL WORK/CHECKING
% ..

ivec = region2imagevec(cl);

clout = region(ivec, 'contiguous_regions', 'noverbose');

% ----------------------------------------------
% Add other fields
% ----------------------------------------------

clindx = ivec.volInfo.cluster;  % wh are in new clusters
nvox = length(clindx);

dat = cat(2, cl.all_data);

if size(dat, 2) == nvox
    
    % matching field; get dat for new clout
    for j = 1:length(clout)
        clout(j).all_data = dat(:, clindx == j);
        
        clout(j).dat = nanmean(clout(j).all_data', 1)';
    end
    
elseif isempty(dat)
    % do nothing
else
    disp('Warning!  all_dat field is wrong size.');
    
end

end



