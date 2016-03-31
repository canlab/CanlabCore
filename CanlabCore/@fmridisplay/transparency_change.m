function transparency_change(o2, multval)
% Change the transparency of blobs in an fmridisplay object
%
% :Inputs:
%
%   **multval:**
%        multiply transparency values by this.
%
% values < 1 makes blobs more transparent, > 1 makes blobs more opaque

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
