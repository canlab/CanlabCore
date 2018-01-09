function cl = merge(cl, wh_merge)
% Merge two or more regions together in a region object.
% Combines fields from all clusters in the named series with the first one
% in the series.
%
% :Usage:
% ::
%
%    wh_merge = [3 4];
%    cl = merge(cl, wh_merge)
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



