function tbl=label_table(atlas_obj)
% label_table Create a table from an atlas object's label fields.
%
% Take an atlas object and produce a MATLAB table whose rows correspond
% to atlas regions and whose columns are label_descriptions, labels,
% labels_2, labels_3, labels_4, and labels_5. Empty auxiliary label
% fields are padded with NaN so that all columns have the same length.
%
% :Usage:
% ::
%
%     tbl = label_table(atlas_obj)
%
% :Inputs:
%
%   **atlas_obj:**
%        An atlas-class object with the following properties:
%
%        - label_descriptions: Descriptions of the labels
%        - labels: Primary labels
%        - labels_2: Secondary labels
%        - labels_3: Tertiary labels
%        - labels_4: Quaternary labels
%        - labels_5: Quinary labels
%
% :Outputs:
%
%   **tbl:**
%        MATLAB table with columns label_descriptions, labels, labels_2,
%        labels_3, labels_4, and labels_5.
%
% :Examples:
% ::
%
%     atlas_obj = load_atlas('canlab2024');
%     tbl = label_table(atlas_obj);
%     disp(tbl);
%
% :See also:
%   - load_atlas
%   - select_atlas_subset
%
% ..
%    Michael Sun, Ph.D. 06/03/2024
% ..

    % Check if the input is valid
    if ~strcmp(class(atlas_obj), 'atlas')
        error('label_table:InvalidInput', 'Input must be an atlas.');
    end

    % Check if any of the labels are empty and pad with NaN
    fields_to_check = {'labels_2', 'labels_3', 'labels_4', 'labels_5'};
    for i = 1:length(fields_to_check)
        if ~isempty(atlas_obj.(fields_to_check{i}))
            atlas_obj.(fields_to_check{i}) = padwithnan(atlas_obj.(fields_to_check{i}), atlas_obj.labels, 2);
        else
            atlas_obj.(fields_to_check{i}) = repmat({''}, 1, numel(atlas_obj.labels));
        end
    end

    tbl=array2table([atlas_obj.label_descriptions, atlas_obj.labels', atlas_obj.labels_2', atlas_obj.labels_3', atlas_obj.labels_4', atlas_obj.labels_5'], ...
        'VariableNames',{'label_descriptions', 'labels', 'labels_2', 'labels_3', 'labels_4', 'labels_5'});


end