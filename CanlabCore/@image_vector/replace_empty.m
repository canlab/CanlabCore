function obj = replace_empty(obj, varargin)
% Replace empty/missing values in an image data object
%
% :Usage:
% ::
%
%    obj = replace_empty(obj, [optional keywords])
%
% Replace missing values in obj.dat stored in obj.removed_voxels and
% obj.removed_images with zeros.  This returns obj.dat in a format that can
% be reconstructed into a 3-D or 4-D image matrix for brain visualization.
%
% :Optional keywords:
%
%   **'voxels' or 'images':**
%        replace only missing voxels/images
%
% ..
%    Tor Wager, 12/1/10
% ..
%
% :See also: remove_empty, zeroinsert, nanremove, naninsert

dovoxels = 1;
doimages = 1;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'voxels', doimages = 0;
            case 'images', dovoxels = 0;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

len = length(obj.removed_images);

% * need proper length checking
% if len > 1 && len ~= size(obj.dat, 2)
%     error('Wrong length for removed_images field!');
% end

    
if doimages && any(obj.removed_images) 
    
    obj.dat = zeroinsert(obj.removed_images, obj.dat')';
    obj.removed_images = false(size(obj.removed_images));

end

if dovoxels && any(obj.removed_voxels)
    
    obj.dat = zeroinsert(obj.removed_voxels, obj.dat);
     
    if isa(obj, 'statistic_image')
%         if ~isempty(obj.p), obj.p = zeroinsert(obj.removed_voxels, obj.p); end % Wani deleted this line, Aug 2012
        if ~isempty(obj.p), obj.p = oneinsert(obj.removed_voxels, obj.p); end  % Wani added this line, Aug 2012 - because the threshold.m function recognizes zeros for p values as significant voxels. 
        if ~isempty(obj.ste), obj.ste = zeroinsert(obj.removed_voxels, obj.ste); end
        if ~isempty(obj.sig), obj.sig = zeroinsert(obj.removed_voxels, obj.sig); end
        
        % 'N' field can be a scalar, indicating number of observations
        % (images), or can be a obs by vox matrix, indicating number of
        % observations at each each voxel. Only do the "remove empty" operation
        % if 'N' is a matrix; don't do it if N is a scalar.
        if ~isscalar(obj.N)
            if ~isempty(obj.N), obj.N = zeroinsert(obj.removed_voxels, obj.N); end
        end
    end
    
    if isa(obj, 'atlas')
        
        has_pmaps = ~isempty(obj.probability_maps) && size(obj.probability_maps, 2) == num_regions(obj);
        
        if has_pmaps
            
            obj.probability_maps = zeroinsert(obj.removed_voxels, obj.probability_maps);
            
        end
        
    end
    
    obj.removed_voxels = false(size(obj.removed_voxels));
    
end


end
