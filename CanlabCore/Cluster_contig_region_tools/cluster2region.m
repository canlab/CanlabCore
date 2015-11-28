function obj = cluster2region(cl)
% Transform a CANlab/SCANlab "clusters" structure into a region object, the
% standard in 2011 toolbox functions and beyond.
%
% ..
%    Tor Wager, Feb 2011
% ..


N = fieldnames(cl);

obj = region;
No = fieldnames(obj);
obj(1:length(cl)) = obj(1);

for k = 1:length(cl)
    
    for i = 1:length(N)
        
        % Look for a field (attribute) with the input name
        wh = strmatch(N{i}, No, 'exact');
        
        if ~isempty(wh)
            
            obj(k).(N{i}) = cl(k).(N{i});
            
        else
            % make mapping with named fields in obj
            
            switch N{i}
                
                case {'from_label', 'from_cluster'}
                    
                    obj(k).custom_info1 = cl(k).(N{i});
                    obj(k).custom_info1_descrip = N{i};
                    
                case 'Z'
                    obj(k).val = cl(k).Z;
                    obj(k).Z = cl(k).Z;
                    
                case {'P', 'image_names'}
                    obj(k).source_images{1} = cl(k).(N{i});
                    
                case {'average_data' 'timeseries' 'contrast_data'}
                    obj(k).dat = cl(k).(N{i});
                    
                case {'threshold'}
                    % do nothing
                    
                otherwise
                    if k == 1
                        warning('cl2region:no comparable field %s',  N{i});
                    end
            end
            
            
        end % case field
        
    end
    
end % cl

end % function



