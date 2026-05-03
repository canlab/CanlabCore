function atlas_obj = atlas_add_L_R_to_labels(atlas_obj)
% atlas_add_L_R_to_labels Standardize lateralization suffixes on atlas labels.
%
% Remove some strings indicating lateralization from atlas labels and add
% new _L and _R suffixes for lateralized regions. The strings stripped from
% the labels are: L_ R_ Left_ Right_ _L _R _Left _Right. New suffixes are
% determined by the modal sign of the x-coordinate for voxels in each
% region (negative = left, positive = right). Regions whose proportional
% L/R asymmetry is below 0.5 are left without a suffix.
%
% :Usage:
% ::
%
%     atlas_obj = atlas_add_L_R_to_labels(atlas_obj)
%
% :Inputs:
%
%   **atlas_obj:**
%        An atlas-class object whose .labels field will be modified.
%
% :Outputs:
%
%   **atlas_obj:**
%        Atlas object with relabeled .labels (with _L/_R suffixes added
%        as appropriate).
%
% :See also:
%   - split_atlas_by_hemisphere
%   - atlas2region

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
