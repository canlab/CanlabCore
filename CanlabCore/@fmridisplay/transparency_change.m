function transparency_change(o2, multval)
% transparency_change Scale the transparency of all blobs in an fmridisplay object.
%
% Multiplies the AlphaData of every blob graphics handle in
% o2.activation_maps by multval, clipping any resulting values above 1.
% Values < 1 make blobs more transparent, > 1 make blobs more opaque.
%
% :Usage:
% ::
%
%     transparency_change(o2, multval)
%
% :Inputs:
%
%   **o2:**
%        An fmridisplay object with one or more activation_maps
%        containing blobhandles.
%
%   **multval:**
%        Scalar to multiply transparency values by. Values < 1 make
%        blobs more transparent, > 1 make blobs more opaque.
%
% :Outputs:
%
%   None. Blob graphics objects in o2 are modified in place.
%
% :See also:
%   - fmridisplay
%   - addblobs
%   - removeblobs

for m = 1:length(o2.activation_maps)
    
    h = o2.activation_maps{m}.blobhandles;
    
    for i = 1:length(h)
        
        d1 = get(h(i), 'AlphaData');
        
        d1 = d1 .* multval;
        
        d1(d1 > 1) = 1;
        
        set(h(i), 'AlphaData', d1);
        
    end
    
end

end
