% I'll have to put these all together and graph them
% for HRF.region

% for HRF.types

function avgHRF = HRF_avg(HRF_cell_array, at)
    % Check if all HRFs have the same atlas

    % Suppress all warnings for now, since they'll mostly flood the screen
    % with Invalid Heights from findpeaks()
    warning('off', 'all');

    % Get the .atlas value of the first structure (assuming HRF_cell_array is not empty)
    first_atlas = HRF_cell_array{1}.atlas;
    first_region = HRF_cell_array{1}.region;
    first_types = HRF_cell_array{1}.types;
    first_conds = HRF_cell_array{1}.params.CondNames;
    
    % Check if all structures have the same .atlas value
    all_same_atlas = true;
    all_same_regions = true;
    all_same_types = true;
    all_same_conds = true;
    
    for i = 2:numel(HRF_cell_array)
        if ~isequal(HRF_cell_array{i}.atlas, first_atlas)
            all_same_atlas = false;
        end
        if ~isequal(HRF_cell_array{i}.region, first_region)
            all_same_regions = false;
        end

        if ~isequal(HRF_cell_array{i}.types, first_types)
            all_same_types = false;
        end

        if ~isequal(HRF_cell_array{i}.params.CondNames, first_conds)
            all_same_conds = false;
        end

    end

    if all_same_atlas
        avgHRF.atlas=first_atlas;
    else
        error('HRF Structures have different .atlas values.\n');
    end

    if all_same_regions
        avgHRF.region=first_region;
    else
        error('HRF Structures have different .region values.\n');
    end

    if all_same_types
        avgHRF.types=first_types;
    else
        error('HRF Structures have different .types values.\n');
    end

    if all_same_conds
        avgHRF.params.CondNames=first_conds;
    else
        error('HRF Structures have different .params.CondNames values.\n');
    end

    % Extract the number of fits, regions, and conditions
    numFits = numel(avgHRF.types);
    numRegions = numel(avgHRF.region);
    numConditions = numel(avgHRF.params.CondNames);
    
    % Loop through fits, regions, and conditions
    for fitIndex = 1:numFits
        
        for regionIndex = 1:numRegions
            
            for conditionIndex = 1:numConditions
                conditionName = avgHRF.params.CondNames{conditionIndex};
                
                % Initialize an empty array to store 'model' vectors
                models = [];
                numPeaks = 0; % Initialize the count of peak elements
                numTroughs = 0; % Initialize the count of peak elements
                
                % Iterate through the structures and extract 'model' vectors
                for dataIndex = 1:numel(HRF_cell_array)
                    modelVector = HRF_cell_array{dataIndex}.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).model;
                    models = [models; modelVector]; % Concatenate the vectors

                    peakVector = HRF_cell_array{dataIndex}.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).peaks;
                    numPeaks = numPeaks + numel(peakVector); % Count the peak elements
                    troughVector = HRF_cell_array{dataIndex}.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).troughs;
                    numTroughs = numTroughs + numel(troughVector); % Count the trough elements

                end
                
                % Compute the average vector for the current condition, region, and fit
                avgVector = mean(models, 1);
                % Compute the average number of elements in the 'peak' field
                avgNumPeaks = numPeaks / numel(HRF_cell_array);
                avgNumTroughs = numTroughs / numel(HRF_cell_array);


                % Store the average vector in the avgStruct
                avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).model = avgVector;
                avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).avgpeaks=avgNumPeaks;
                avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).avgtroughs=avgNumTroughs;

                % Store the constituent vectors in the avgStruct
                avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).models = models;

                [avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).peaks, avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).troughs]=detectPeaksTroughs(avgVector', false);
                [~, region_numVox, ~, ~]=at.select_atlas_subset(avgHRF.region(regionIndex), 'exact').get_region_volumes;
                [avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).peaks_voxnormed, avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).troughs_voxnormed]=detectPeaksTroughs(avgVector'/region_numVox, false);

                % Number of phases
                start_times = [avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).peaks.start_time, avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).troughs.start_time];
                end_times = [avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).peaks.end_time, avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).troughs.end_time];
                phases = [start_times(:), end_times(:)];
                unique_phases = unique(phases, 'rows');
                avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phases = mat2cell(unique_phases, ones(size(unique_phases, 1), 1), 2);
    
                % Loop over phases
                for p = 1:numel(avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phases)
                    features = {'peaks', 'troughs'};
                    for f = 1:2
                        feat = features{f};
                        
                        % Count the number of features (peaks or troughs)
                        start_times = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).(feat).start_time});
                        current_phase_start = avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phases{p}(1);
                        avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phase(p).(feat) = sum(start_times == current_phase_start);
                        
                        if avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phase(p).(feat) > 0
                            % display(['Estimating ' , num2str(fitIndex), '_', rois{r}, '_', first_conds{conditionIndex}, '_', ' Now...!'])
                            % display(['Phase ' , num2str(p), 'Feature ', feat])
                            idx = find(start_times == current_phase_start);
                            auc = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).(feat).AUC});
                            height = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).(feat).height});
                            time_to_peak = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).(feat).time_to_peak});
                            half_height = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).(feat).half_height});
                            
                            feat_voxnormed = strcat(feat, '_voxnormed');
                            auc_voxnormed = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).(feat_voxnormed).AUC});
                            height_voxnormed = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).(feat_voxnormed).height});
                            time_to_peak_voxnormed = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).(feat_voxnormed).time_to_peak});
                            half_height_voxnormed = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).(feat_voxnormed).half_height});
                            
                            avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phase(p).auc = unique(auc(idx));
                            avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phase(p).auc_voxnormed = unique(auc_voxnormed(idx));
    
                            avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phase(p).height = height(idx);
                            avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phase(p).height_voxnormed = height_voxnormed(idx);
                            avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phase(p).time_to_peak = time_to_peak(idx);
                            avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phase(p).time_to_peak_voxnormed = time_to_peak_voxnormed(idx);
                            avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phase(p).half_height = half_height(idx);
                            avgHRF.fit{fitIndex}{regionIndex}.(first_conds{conditionIndex}).phase(p).half_height_voxnormed = half_height_voxnormed(idx);
                        end
                    end
                end
            end
        end
    end
end
    
% 
% HRF{file}.fit{1}{1}.hot.model
% HRF{file}.fit{1}{1}.warm.model
% HRF{file}.fit{1}{1}.imagine.model
% 
% 
% 
% HRF{file}.fit{1}{1}.model.h


