
% RUn [atlas_of_subdivisions, region_values] = subdivide_by_atlas(activation_map, atl, [optional inputs])

dosep = true;
donames = true;


% Separate subregions with positive and negative values if requested
% -------------------------------------------------------------------------
if dosep
    % separate pos and neg
    [poscl, negcl] = posneg_separate(cl);
    
    cl = [poscl negcl];
    ispos = [true(1, length(poscl)) false(1, length(negcl))]; % logical for splitting combined cl later
    
    clear poscl negcl
    
    fprintf('\n%s\nPositive Effects\n', sep_str)
else
    %     % just return cl in poscl
    %     poscl = cl;
    %     negcl = [];
    fprintf('\n%s\nTable of all regions\n', sep_str)
end



    if donames
        % we have already entered names
        Region = table(region_table.Region, 'VariableNames', {'Region'});
    else
        % use modal label
        Region = table(region_table.modal_label, 'VariableNames', {'Region'});
    end
    
    Volume = table(region_table.Region_Vol_mm, 'VariableNames', {'Volume'});
    Atlas_coverage = region_table(:, [6 7 4]);
    XYZ = table(round(cat(1, cl.mm_center)), 'VariableNames', {'XYZ'});
    
    % Z = get_max_Z(cl);
    Z = get_signed_max(cl, 'Z', 'maxZ');  % use function because may be empty, handle if so
    
    results_table = [Region Volume XYZ Z Atlas_coverage];
    results_table.region_index = (1:size(region_table, 1))';
    
    results_table_pos = results_table(ispos, :);
    results_table_neg = results_table(~ispos, :);
    
    % Sort, if asked for (default = yes)
    if dosortrows
    
        % Replace empty strings so sort will work
        whempty = cellfun(@isempty, results_table_pos.modal_label_descriptions);
        results_table_pos.modal_label_descriptions(whempty) = {'X_No_label'};
        
        whempty = cellfun(@isempty, results_table_neg.modal_label_descriptions);
        results_table_neg.modal_label_descriptions(whempty) = {'X_No_label'};
        
        results_table_pos = sortrows(results_table_pos, 'modal_label_descriptions');
        results_table_neg = sortrows(results_table_neg, 'modal_label_descriptions');
        
        % Manual - not as good because Matlab's table methods handle this well.
        %             % Replace empty and get unique labels to sort by
        %             whempty = cellfun(@isempty, results_table_pos.modal_label_descriptions);
        %             results_table_pos.modal_label_descriptions(whempty) = {'No_label'};
        %             u = unique(results_table_pos.modal_label_descriptions);
        % 
        %             [~, ~, condf] = string2indicator(results_table_pos.modal_label_descriptions)
        % 

    end
    
    
    % Now split into positive and neg sub-tables and display
    
    if any(ispos)
        disp(results_table_pos)
    else
        disp('No regions to display');
    end
    
    fprintf('\nNegative Effects\n')
    if any(~ispos)
        disp(results_table_neg)
    else
        disp('No regions to display');
    end
    
end

if dolegend == false || isempty(table_legend_text)
    return
end

% clean up text
if length(table_legend_text) > 2
    table_legend_text(2:3) = [];
end