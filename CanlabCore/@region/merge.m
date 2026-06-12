function cl = merge(cl, wh_merge)
% merge Merge two or more regions together within a region object.
%
% Combine the regions indexed by wh_merge into a single region (the
% first one in wh_merge), updating fields as appropriate: text fields
% are concatenated with ' MERGED WITH ' separators, voxel-list fields
% (XYZ, XYZmm, Z, threshold, all_data) are concatenated horizontally,
% the .val field is concatenated vertically, and per-region averages
% (timeseries, contrastdata, dat) are weighted-averaged by voxel count.
% Center and mm_center are recomputed and the merged regions other
% than the first are removed from the array.
%
% :Usage:
% ::
%
%     wh_merge = [3 4];
%     cl = merge(cl, wh_merge)
%
% :Inputs:
%
%   **cl:**
%        A region-class object array.
%
%   **wh_merge:**
%        Vector of indices into cl identifying the regions to merge.
%        The first entry receives the merged result; the remaining
%        entries are deleted from cl.
%
% :Outputs:
%
%   **cl:**
%        Region object array with the merged region in position
%        wh_merge(1) and the other merged elements removed.
%
% :See also:
%   - region
%   - reparse_continguous
%
% ..
%    Tor Wager, April 2011
% ..

N = fieldnames(cl(wh_merge(1))); % add stuff to the first one in the list

for j = 2:length(wh_merge)
    
    for i = 1:length(N)
        
        try
            
            switch N{i}
                
                % Update as strings by concatenating
                case {'title' 'shorttitle' 'descrip1' 'descrip2' 'custom_info1' 'custom_info2' 'custom_info1_descrip' 'custom_info2_descrip'}
                    
                    cl(wh_merge(1)).(N{i}) = [cl(wh_merge(1)).(N{i}) ' MERGED WITH ' cl(wh_merge(j)).(N{i})];
                    
                    % Merge by concatenating vertically
                case {'val'}
                    cl(wh_merge(1)).(N{i}) = [cl(wh_merge(1)).(N{i}); cl(wh_merge(j)).(N{i})];
                    
                    % Merge by concatenating horizontally
                case {'XYZ' 'XYZmm' 'Z' 'threshold' 'all_data'}
                    cl(wh_merge(1)).(N{i}) = [cl(wh_merge(1)).(N{i}) cl(wh_merge(j)).(N{i})];
                    
                    % Merge  by weighted average of columns/data
                case {'timeseries' 'contrastdata' 'dat'}
                    
                    mysum = cl(wh_merge(1)).numVox .* cl(wh_merge(1)).(N{i}) + cl(wh_merge(j)).numVox .* cl(wh_merge(j)).(N{i});
                    mysum = mysum ./ (cl(wh_merge(1)).numVox + cl(wh_merge(j)).numVox);
                    
                    cl(wh_merge(1)).(N{i}) = mysum;
                    
                    % Don't update these now
                case {'val_descrip' 'numpeaks' 'center' 'mm_center' 'source_images'}
                    
                    
                    % Check for compatibility
                case {'M' 'dim' 'voxSize'}
                    
                    mycomp = cl(wh_merge(1)).(N{i}) - cl(wh_merge(j)).(N{i});
                    
                    if any(mycomp(:))
                        fprintf('WARNING! Field %s does not match for regions %3.0f and %3.0f.  Combined values may be meaningless.\n', N{i}, 1, j);
                    end
                    
                    
            end % case
            
        catch
            fprintf('WARNING! Field %s has sizes that do not match for regions %3.0f and %3.0f.  Not combining this field.\n', N{i}, 1, j);
            
        end
        
    end % field i
    
    
    % Update center, mm_center
    cl(wh_merge(1)).center = nanmean(cl(wh_merge(1)).XYZ, 2)';
    cl(wh_merge(1)).mm_center = nanmean(cl(wh_merge(1)).XYZmm, 2)';
    
end % cluster j

cl(wh_merge(1)).numVox = sum(vertcat(cl(wh_merge).numVox));
cl(wh_merge(2:end)) = [];

end % function



