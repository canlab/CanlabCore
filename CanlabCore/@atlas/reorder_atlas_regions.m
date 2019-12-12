function obj = reorder_atlas_regions(obj, wh_order)
% Reorder a set of regions in an atlas object
%
% Reorder a set of regions in an atlas object according to an integer
% vector, wh_order, that specifies the original region numbers that should
% become regions 1...n.  e.g., if wh_order is [3 1], region 3 becomes
% region 1, region 1 becomes region 2.
%
% obj = reorder_atlas_regions(obj, wh_order)

% Programmers' notes:
% Created by Tor Wager
% 12/2019: Tor fixed bug in reordering; probability_maps field was not
% preivously being considered, which resulted in wrong ordering when
% present. .dat is rebuilt from probability_maps if available.

n = num_regions(obj);

if any(wh_order > n)
    error('wh_order contains integer values that do not have corresponding region indices.');
end


newdat = int32(zeros(size(obj.dat)));

for i = 1:length(wh_order)
    
    newdat(obj.dat == wh_order(i)) = i;
    
end

obj.dat = newdat;

% Reorder probability maps
if ~isempty(obj.probability_maps)
    
    if ~size(obj.probability_maps, 2) == n
        error('obj.probabilty_maps is not empty, but is the wrong size! illegal object.');
    end
    
    obj.probability_maps = obj.probability_maps(:, wh_order);
    
end

% Reorder labels
myfields = {'labels' 'labels_2' 'labels_3' 'labels_4' 'labels_5' 'property_descriptions', 'label_descriptions'};

for i = 1:length(myfields)
    
    if ~isempty(obj.(myfields{i})) && length(obj.(myfields{i})) == n
        
        obj.(myfields{i}) = obj.(myfields{i})(wh_order);
        
    end
    
end

end % function

