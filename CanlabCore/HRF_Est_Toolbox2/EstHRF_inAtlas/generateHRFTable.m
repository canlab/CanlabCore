function T=generateHRFTable(HRF)
    % Generates a summary table from an HRF Structure.
    % Michael Sun, Ph.D.
    % - Takes HRF Structure object generated from EstimateHRF_inAtlas()
    %
    % *Usage:
    % ::
    %    T=generateHRFTable(HRF_structure)

    % Pre-allocate a cell array to store the table data
    tblData = cell(0, 10); % You may need to adjust the number of columns based on the data you're storing
    
    % Iterate through the nested loops as before
    for c = 1:numel(HRF.CondNames)
        for r = 1:numel(HRF.region)
            for t = 1:numel(HRF.types)
            
                %disp(['Number of phases: ', num2str(numel(HRF.fit{t}{r}.(HRF.CondNames{c}).phases))]);
                %for p = 1:numel(HRF.fit{t}{r}.(HRF.CondNames{c}).phases)
                    
                % Define the basic information
                type = HRF.types{t};
                region = HRF.region{r};
                condition = HRF.CondNames{c};
                phase_num = 0;
                phase_span = [0 0];
                peaks = 0;
                troughs = 0;
                time_to_peak = 0;
                auc = 0;
                auc_voxnormed = 0;
    
    
                if numel(HRF.fit{t}{r}.(HRF.CondNames{c}).phases) > 1
                    for p = 1:numel(HRF.fit{t}{r}.(HRF.CondNames{c}).phases)
            
                        phase_num = p;
                        phase_span = HRF.fit{t}{r}.(HRF.CondNames{c}).phases{p};
                        
                        % Check for peaks and troughs
                        peaks = HRF.fit{t}{r}.(HRF.CondNames{c}).phase(p).peaks;
                        troughs = HRF.fit{t}{r}.(HRF.CondNames{c}).phase(p).troughs;
                        %disp(['Peaks: ', num2str(peaks), ', Troughs: ', num2str(troughs)]);
                        % Extract other data
                        time_to_peak = HRF.fit{t}{r}.(HRF.CondNames{c}).phase(p).time_to_peak;
                        auc = HRF.fit{t}{r}.(HRF.CondNames{c}).phase(p).auc;
                        auc_voxnormed = HRF.fit{t}{r}.(HRF.CondNames{c}).phase(p).auc_voxnormed;
                        % Append the data to the cell array
                        tblData = [tblData; {type, region, condition, phase_num, phase_span, peaks, troughs, time_to_peak, auc, auc_voxnormed}];
                    end
                else
                    % Append the data to the cell array
                    %disp(['c: ', num2str(c), ', t: ', num2str(t), ', r: ', num2str(r)]);
                    tblData = [tblData; {type, region, condition, phase_num, phase_span, peaks, troughs, time_to_peak, auc, auc_voxnormed}];
    
                end
            end
        end
    end
    
    % Convert the cell array to a table
    T = cell2table(tblData, 'VariableNames', {'Type', 'Region', 'Condition', 'PhaseNum', 'PhaseSpan', 'Peaks', 'Troughs', 'TimeToPeak', 'AUC', 'AUC_VoxNormed'});
    
    % Display the table
    disp(T);
end