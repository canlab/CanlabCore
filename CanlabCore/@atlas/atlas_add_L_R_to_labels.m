function atlas_obj = atlas_add_L_R_to_labels(atlas_obj)
% Removes some strings indicating lateralization from atlas labels and adds new _L and _R suffixes for lateralized regions.
% - strings replaced are L_ R_ Left_ Right_ _L _R _Left _Right

labels = atlas_obj.labels;

mypatterns = {'L_' 'R_' 'Left_' 'Right_' '_L' '_R' '_Left' '_Right'};

for i = 1:length(mypatterns)
    labels = regexprep(labels, mypatterns{i}, '');
end

mysuffix = get_suffix(atlas_obj);

for i = 1:length(labels)
    
    labels{i} = [labels{i} mysuffix{i}];
    
end

atlas_obj.labels = labels;

end % function






function mysuffix = get_suffix(atlas_obj)

r = atlas2region(atlas_obj);

k = length(r);

for j = 1:k
    
    voxsign = sign(r(j).XYZmm(1, :)); % sign of each x coord. - = left
    modal_x = mode(voxsign);
    
    switch modal_x
        case -1
            labelstr = '_L';
        case 1
            labelstr = '_R';
        case 0
            labelstr = '';
    end
    
    % proportional asymmetry. 1 is vary lateralized. 0 is balanced
    % across L/R.  If symmetrical, then label _M
    prop_asym = abs(sum(voxsign == -1) - sum(voxsign == 1)) ./ length(voxsign);
    
    if prop_asym < .5, labelstr = ''; end
    
    mysuffix{j} = labelstr;
    
    
end % j

end % subfunction
