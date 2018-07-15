function [obj, has_pmaps, has_index, missing_regions] = check_properties(obj, varargin)
% Check properties and enforce some variable types
%
% obj = check_properties(obj)
%
% Optional arguments:
% 'compress_index' : if index numbers are not consecutive integers, some functions,
%                    like select_atlas_subset and num_regions, will not work
%                    This rebuilds the index.  This could happen after resampling or masking.
%                    If you have probability maps entered, probability_maps_to_region_index will do
%                    this. But if not, you may want to do this.

% Defaults

compress_index = false;

% Optional inputs

if any(strcmp(varargin, 'compress_index')), compress_index = true; end

% Count regions
[n_regions, n_regions_with_data, missing_regions] = num_regions(obj);

% Data
has_pmaps = ~isempty(obj.probability_maps) && size(obj.probability_maps, 2) == n_regions;
has_index = ~isempty(obj.dat);

% enforce sparse double prob maps
if has_pmaps
    
    obj.probability_maps = sparse(double(obj.probability_maps));
    
end

% Create index if needed
if has_pmaps && ~has_index
    obj = probability_maps_to_region_index(obj);
end

% could check for valid data here and fix, warn, or error

valid_pmaps = all(obj.probability_maps(:) >= 0) && all(obj.probability_maps(:) <= 1);
valid_index = size(obj.dat, 2) == 1 && all(abs(obj.dat(:) - round(obj.dat(:))) < 100 * eps);


% Check integers
u = unique(double(obj.dat(obj.dat ~= 0)));

if ~all(u == round(u))
    warning('Some region values on obj.dat are not integers. Adjusting.');
    obj.dat = round(obj.dat);
end

% n_regions_from_indx = length(u);


% Enforce variable types
% ----------------------------------------------------------------------------

obj.probability_maps = sparse(double(obj.probability_maps));
obj.dat = int32(full(round(obj.dat)));

myfields = {'labels' 'labels_2' 'labels_3' 'labels_4' 'labels_5' 'property_descriptions', 'label_descriptions'};

for i = 1:length(myfields)
    
    if ~iscell(obj.(myfields{i})), obj.(myfields{i}) = {obj.(myfields{i})}; end
    
end


% Check for names/labels and create placeholders if missing
% ----------------------------------------------------------------------------

if isempty(obj.label_descriptions)
    obj.label_descriptions = cell(n_regions, 1);
end

if isempty(obj.labels)
    obj.labels = cell(1, n_regions);
end

for i = 1:n_regions
    
    if length(obj.labels) < i || isempty(obj.labels{i})
        
        obj.labels{i} = sprintf('Reg_%d', i);
        
    end
    
    if length(obj.label_descriptions) < i || isempty(obj.label_descriptions{i})
        
        obj.label_descriptions{i} = sprintf('Region %d from %s', i, obj.atlas_name);
        
    end
    
end

% Enforce row/column format (avoids compatibility problems later)

if iscolumn(obj.labels), obj.labels = obj.labels'; end

if ~iscolumn(obj.label_descriptions), obj.label_descriptions = obj.label_descriptions'; end



% Deal with missing integers
% -------------------------------------------------------------------------
% Ok to have some missing integers. Can happen if resampling spaces.
% but some functions, like select_atlas_subset and num_regions, will not work if index
% values are off!
% if ~(max(obj.dat) == n_regions) || n_regions_from_indx ~= n_regions
%     warning('Integer index in obj.dat should have one integer value per label, without skipping any integers!!');
% end

if n_regions ~= n_regions_with_data || any(missing_regions)
    
    if compress_index
        % disp('Fixing');
        
        newdat = obj.dat;
        
        for i = 1:length(u)
            
            newdat(obj.dat == u(i)) = i;
            
        end
        
        obj.dat = newdat;
        
        if n_regions == n_regions_with_data
            % We have missing regions but a complete set of labels
            % do not adjust labels
            return
        end
        
        % Trim labels
        
        for i = 1:length(myfields)
            
            if length(obj.(myfields{i})) == n_regions  % labels exist but need trimming
                obj.(myfields{i}) = obj.(myfields{i})(u);
            end
            
        end
        
    else
        % Just print warning, we are not compressing
        
        warning('Integer index in obj.dat should have one integer value per label, without skipping any integers!!');
        disp(missing_regions')
        disp('Run check_properties with ''compress_index'' option to fix.');
        
    end
    
end


end % function

