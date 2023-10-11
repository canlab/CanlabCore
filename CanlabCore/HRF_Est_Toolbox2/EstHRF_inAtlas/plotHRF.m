function plotHRF(HRF, t, r)
    % Generates a plot from an HRF Structure given a specified fit type and a region name.
    % Michael Sun, Ph.D.
    % - Takes HRF Structure object generated from EstimateHRF_inAtlas()
    % - t is a cellstring for fit type e.g., 'FIR', 'IL', or 'CHRF'
    % - r is a cellstring for region label from atlas.labels. e.g., 'CA2_Hippocampus_' 
    %
    % *Usage:
    % ::
    %    plotHRF(HRF_structure, 'FIR', 'CA2_Hippocampus_')

    % Check if the specified region 'r' and type 't' exist in HRF
    if ~ismember(r, HRF.region)
        error('Invalid region specified.');
    end
    
    if ~ismember(t, HRF.types)
        error('Invalid type specified.');
    end

    % Find the indices of the specified region and type
    reg = find(strcmp(HRF.region, r));
    typ = find(strcmp(HRF.types, t));

    % Get the list of conditions from HRF_PARAMS
    conds = HRF.params.CondNames;

    % Initialize a cell array to store the model matrices
    array3D = cell(1, length(conds));

    % Populate the cell array with model matrices for each condition
    for f = 1:numel(conds)
        if isfield(HRF.fit{typ}{reg}, conds{f}) && isfield(HRF.fit{typ}{reg}.(conds{f}), 'model')
            array3D{f} = HRF.fit{typ}{reg}.(conds{f}).model;
        else
            error('Model data not found for condition: %s', conds{f});
        end
    end

    % Check if any model matrices were found
    if all(cellfun(@isempty, array3D))
        error('No model data found for the specified type and region.');
    end

    % Concatenate the model matrices along the third dimension
    array3D = cat(3, array3D{:});

    % Initialize color map and legend entries
    % colors = jet(numel(conds));

    % Create a figure
    figure;

    % Loop through each condition and plot
    for cond = 1:numel(conds)
        subplot(numel(conds), 1, cond);

        if isfield(HRF.fit{typ}{reg}, conds{cond}) && isfield(HRF.fit{typ}{reg}.(conds{cond}), 'models')
            for m = 1:height(HRF.fit{typ}{reg}.(conds{cond}).models)
                plot(HRF.fit{typ}{reg}.(conds{cond}).models(m,:), 'Color', [0.7,0.7,0.7], 'LineWidth', 0.2);
                hold on;
            end

        detectPeaksTroughs(squeeze(array3D(:, :, cond))', true);
        hold on;
        hline(0);
        region=format_strings_for_legend(r);
        region=region{1};
        title({['Condition ', conds{cond}], ['Fit-type: ', t], ['Region: ', region]}, 'Interpreter', 'none');
    end


end
