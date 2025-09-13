function tbl=label_table(atlas_obj)
    % LABEL_TABLE Create a table from an atlas object
    %
    %   TBL = LABEL_TABLE(ATLAS_OBJ) takes an atlas object ATLAS_OBJ and
    %   creates a table TBL containing the label descriptions and various
    %   labels.
    %
    %   Input:
    %     - atlas_obj: An object containing atlas data with the following
    %       properties:
    %         * label_descriptions: Descriptions of the labels
    %         * labels: Primary labels
    %         * labels_2: Secondary labels
    %         * labels_3: Tertiary labels
    %         * labels_4: Quaternary labels
    %         * labels_5: Quinary labels
    %
    %   Output:
    %     - tbl: A table containing the label descriptions and labels from
    %       the atlas object, with appropriate variable names.
    %
    % Example:
    %   atlas_obj = load_atlas('canlab2024');
    %   tbl = label_table(atlas_obj);
    %   disp(tbl);
    %
    % Michael Sun, Ph.D. 06/03/2024

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