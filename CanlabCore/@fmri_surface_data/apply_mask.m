function obj = apply_mask(obj, mask, varargin)
% apply_mask Mask a grayordinate object by zeroing out-of-mask grayordinates.
%
% :Usage:
% ::
%     obj = apply_mask(obj, mask)
%     obj = apply_mask(obj, mask, 'invert')
%
% Surface analogue of fmri_data.apply_mask, greatly simplified per design
% decision D5b: because all objects share a standardized grayordinate space,
% masking is intrinsic -- there is no fmri_mask_image, no resampling, and no
% empty-removal. Out-of-mask grayordinate rows are simply set to 0; .dat keeps
% its full size and the geometry is unchanged.
%
% :Inputs:
%   **obj:**  an fmri_surface_data object.
%   **mask:** one of
%       - logical/numeric vector [nGrayordinates x 1] (nonzero = keep)
%       - another fmri_surface_data on the same space (nonzero/non-NaN = keep)
%
% :Optional Inputs:
%   **'invert':** keep the complement (zero the in-mask grayordinates instead).
%
% :Outputs:
%   **obj:** masked object (out-of-mask rows zeroed).
%
% :See also: fmri_surface_data, compare_space

invert = any(strcmpi(varargin, 'invert'));

if isa(mask, 'fmri_surface_data')
    if compare_space(obj, mask) ~= 0
        error('fmri_surface_data:apply_mask:space', ...
            ['Mask is not on the same grayordinate space as the data ' ...
             '(compare_space ~= 0). Resample first.']);
    end
    keep = any(mask.dat ~= 0 & ~isnan(mask.dat), 2);
elseif islogical(mask) || isnumeric(mask)
    keep = logical(mask(:));
    if numel(keep) ~= size(obj.dat, 1)
        error('fmri_surface_data:apply_mask:length', ...
            'Mask vector length (%d) must equal the number of grayordinates (%d).', ...
            numel(keep), size(obj.dat, 1));
    end
else
    error('fmri_surface_data:apply_mask:type', ...
        'mask must be a logical/numeric vector or an fmri_surface_data object.');
end

if invert, keep = ~keep; end

obj.dat(~keep, :) = 0;
obj.history{end+1} = sprintf('apply_mask: kept %d / %d grayordinates (zeroed the rest)', ...
    nnz(keep), numel(keep));
end
