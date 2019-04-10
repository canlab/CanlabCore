function dat = remove_empty(dat, varargin)
% remove vox: logical vector of custom voxels to remove, VOX x 1
%
% remove im: logical vector of custom images to remove, 1 x IMAGES
%
% indices of removed data will be stored in removed_voxels and
% removed_images fields, to preserve ability to later reconstruct into 3D images
%
% :Usage:
% ::
%
%    dat = remove_empty(dat, [logical vector of custom voxels to remove], [logical vector of imgs to remove])
%
% Indicator vectors stored in:
% removed_images
% removed_voxels
%
% :See also: replace_empty

% force logical
dat.removed_images = logical(dat.removed_images);
dat.removed_voxels = logical(dat.removed_voxels);

% Make sure field orientations conform - diff versions of code are
% different
if length(dat.removed_images) > size(dat.removed_images, 1), dat.removed_images = dat.removed_images'; end
if length(dat.removed_voxels) > size(dat.removed_voxels, 1), dat.removed_voxels = dat.removed_voxels'; end
  
% get any new empty images that didn't exist before
empty_images = all(dat.dat == 0 | isnan(dat.dat), 1)';
empty_voxels = all(dat.dat' == 0 | isnan(dat.dat'), 1)';

% add in previously removed voxels/images
if any(dat.removed_images) 
    empty_images = logical(zeroinsert(dat.removed_images, empty_images));
end

if any(dat.removed_voxels)
    empty_voxels = logical(zeroinsert(dat.removed_voxels, empty_voxels));
end

% 2010a : Getting error with automatic sparse vector creation
empty_images = full(empty_images);
empty_voxels = full(empty_voxels);

% Now empties are original length.

v = length(empty_voxels);
k = length(empty_images);

dat.history{end+1} = sprintf('removed %3.0f empty voxels and %3.0f empty images', sum(empty_voxels), sum(empty_images));


% add custom removal specifications
% -----------------------------------------------------------------------------------

if ~isempty(varargin) && any(varargin{1})
    
    remvox = varargin{1};
    if length(remvox) > size(remvox, 1), remvox = remvox'; end
    
    if length(remvox) ~= v
        error('Input to remove_empty for custom voxels must be logical vector of length == total num voxels')
    end

    % new voxels to remove
    empty_voxels = empty_voxels | remvox;
    
    dat.history{end+1} = sprintf('removed %3.0f custom-specified voxels', sum(remvox));

end

if length(varargin) > 1 && any(varargin{2})
    
    remimgs = varargin{2};
    if length(remimgs) > size(remimgs, 1), remimgs = remimgs'; end
    
    if length(remimgs) ~= k
        error('Input to remove_empty for custom images must be logical vector of length == num images')
    end
    
    % new images to remove
    empty_images = empty_images | remimgs;
    
    dat.history{end+1} = sprintf('removed %3.0f custom-specified images', sum(remimgs));

end

% now do two things: update overall removed vox/image list, and remove new
% empties/to-removes
% -----------------------------------------------------------------------------------

prevvox = dat.removed_voxels;
previmgs = dat.removed_images;

% if removed_voxels is [], replace
if isempty(dat.removed_voxels), dat.removed_voxels = false(size(empty_voxels)); end
if isempty(dat.removed_images), dat.removed_voxels = false(size(empty_images)); end

% keep overall list
dat.removed_voxels = dat.removed_voxels | empty_voxels;
dat.removed_images = dat.removed_images | empty_images;

% remove new voxels. make indices current with current .dat field.
if any(prevvox)
    empty_voxels(prevvox) = [];
    empty_images(previmgs) = [];
end

% remove data
% -----------------------------------------------------------------------------------
if ~isempty(dat.dat)
    dat.dat(empty_voxels, :) = [];
    dat.dat(:, empty_images) = [];
else
    return
end

% Remove  voxels from special statistic_image fields

if isa(dat, 'statistic_image')
    
    % 'N' field can be a scalar, indicating number of observations
    % (images), or can be a obs by vox matrix, indicating number of
    % observations at each each voxel. Only do the "remove empty" operation
    % if 'N' is a matrix; don't do it if N is a scalar.
    if isscalar(dat.N)
        myfield = {'p' 'ste' 'sig'};
    else
        myfield = {'p' 'ste' 'sig' 'N'};
    end
    
    for myfields = myfield
        
        % check orientation
        if isrow(dat.(myfields{1})), dat.(myfields{1}) = dat.(myfields{1})'; end
        
        if ~isempty(dat.(myfields{1}))
            dat.(myfields{1})(empty_voxels, :) = [];
        end
        
    end
    
end

% Remove voxels from special atlas fields

if isa(dat, 'atlas')
    
    has_pmaps = ~isempty(dat.probability_maps) && size(dat.probability_maps, 2) == num_regions(dat);
    
    if has_pmaps
        dat.probability_maps(empty_voxels, :) = [];
    end
    
end

% error checking
% -----------------------------------------------------------------------------------
numbad = sum(dat.removed_voxels);
sy = size(dat.dat, 1);
len = length(dat.removed_voxels);

if numbad + sy ~= len
    disp('Illegal dat.removed_voxels vector. length must equal len of original dataset, size(y, 1) + # left-out cases.');
    error('Left out = %3.0f, size(y, 1) = %3.0f, sum = %3.0f, length(wasbad) = %3.0f', numbad, sy, numbad + sy, len);
end
    

end
