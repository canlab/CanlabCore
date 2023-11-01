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

    conditions=avgHRF.params.CondNames;

    % Extract the number of fits, regions, and conditions
    numFits = numel(avgHRF.types);
    numRegions = numel(avgHRF.region);
    numConditions = numel(avgHRF.params.CondNames);
    
    % Loop through fits, regions, and conditions
    for fitIndex = 1:numFits
        
        for regionIndex = 1:numRegions
            avgVector_across_conditions=[];
            models3d=[];
            for conditionIndex = 1:numConditions
                conditionName = conditions{conditionIndex};
                
                % Initialize an empty array to store 'model' vectors
                models = [];
                phases=[];
                peaks=[];
                troughs=[];
                numPhases = 0;
                numPeaks = 0; % Initialize the count of peak elements
                numTroughs = 0; % Initialize the count of peak elements
                
                % Iterate through the structures and extract 'model' vectors
                for dataIndex = 1:numel(HRF_cell_array)
                    modelVector = HRF_cell_array{dataIndex}.fit{fitIndex}{regionIndex}.(conditionName).model;
                    models = [models; modelVector]; % Concatenate the vectors

                    phaseStructs = HRF_cell_array{dataIndex}.fit{fitIndex}{regionIndex}.(conditionName).phases;
                    numPhases = numel(phaseStructs);
                    phases = [phases; numPhases];

                    peakStructs = HRF_cell_array{dataIndex}.fit{fitIndex}{regionIndex}.(conditionName).peaks;
                    numPeaks = numel(peakStructs); % Count the peak elements
                    peaks = [peaks; numPeaks];

                    troughStructs = HRF_cell_array{dataIndex}.fit{fitIndex}{regionIndex}.(conditionName).troughs;
                    numTroughs = numel(troughStructs); % Count the trough elements
                    troughs = [troughs; numTroughs];
                end
                
                % Compute the average vector for the current condition, region, and fit
                avgVector = mean(models, 1);
                avgVector_across_conditions=[avgVector_across_conditions; avgVector];

                % Compute the standard deviation for the current condition, region, and fit
                stdVector = std(models, 0, 1); % second parameter = 0 normalizes by N-1, which is typically desired when estimating the population standard deviation from a sample
                wseVector=barplot_get_within_ste(models);

                avg_per_subject = mean(models, 2);

                avgPhases = mean(phases, 1);
                stdPhases = std(phases, 0, 1);

                avgPeaks = mean(peaks, 1);
                stdPeaks = std(peaks, 0, 1);

                avgTroughs = mean(troughs, 1);
                stdTroughs = std(troughs, 0, 1);

                models3d=cat(3, models3d, models);

                % Store the average vector in the avgStruct
                avgHRF.fit{fitIndex}{regionIndex}.(conditionName).model = avgVector;
                avgHRF.fit{fitIndex}{regionIndex}.(conditionName).stddev = stdVector;
                avgHRF.fit{fitIndex}{regionIndex}.(conditionName).wse = wseVector;

                avgHRF.fit{fitIndex}{regionIndex}.(conditionName).avgpeaks = avgPeaks;
                avgHRF.fit{fitIndex}{regionIndex}.(conditionName).stdpeaks = stdPeaks;

                avgHRF.fit{fitIndex}{regionIndex}.(conditionName).avgtroughs = avgTroughs;
                avgHRF.fit{fitIndex}{regionIndex}.(conditionName).stdtroughs = stdTroughs;

                % Store the constituent vectors in the avgStruct
                avgHRF.fit{fitIndex}{regionIndex}.(conditionName).models = models;

                [avgHRF.fit{fitIndex}{regionIndex}.(conditionName).peaks, avgHRF.fit{fitIndex}{regionIndex}.(conditionName).troughs]=detectPeaksTroughs(avgVector', false);
                [~, region_numVox, ~, ~]=at.select_atlas_subset(avgHRF.region(regionIndex), 'exact').get_region_volumes;
                [avgHRF.fit{fitIndex}{regionIndex}.(conditionName).peaks_voxnormed, avgHRF.fit{fitIndex}{regionIndex}.(conditionName).troughs_voxnormed]=detectPeaksTroughs(avgVector'/region_numVox, false);

                % Number of phases
                start_times = [avgHRF.fit{fitIndex}{regionIndex}.(conditionName).peaks.start_time, avgHRF.fit{fitIndex}{regionIndex}.(conditionName).troughs.start_time];
                end_times = [avgHRF.fit{fitIndex}{regionIndex}.(conditionName).peaks.end_time, avgHRF.fit{fitIndex}{regionIndex}.(conditionName).troughs.end_time];
                phases = [start_times(:), end_times(:)];
                unique_phases = unique(phases, 'rows');
                avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phases = mat2cell(unique_phases, ones(size(unique_phases, 1), 1), 2);
    
                % Loop over phases
                for p = 1:numel(avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phases)
                    features = {'peaks', 'troughs'};
                    for f = 1:2
                        feat = features{f};
                        
                        % Count the number of features (peaks or troughs)
                        start_times = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(conditionName).(feat).start_time});
                        current_phase_start = avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phases{p}(1);
                        avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phase(p).(feat) = sum(start_times == current_phase_start);
                        
                        if avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phase(p).(feat) > 0
                            % display(['Estimating ' , num2str(fitIndex), '_', rois{r}, '_', conditionName, '_', ' Now...!'])
                            % display(['Phase ' , num2str(p), 'Feature ', feat])
                            idx = find(start_times == current_phase_start);
                            auc = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(conditionName).(feat).AUC});
                            height = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(conditionName).(feat).height});
                            time_to_peak = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(conditionName).(feat).time_to_peak});
                            half_height = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(conditionName).(feat).half_height});
                            
                            feat_voxnormed = strcat(feat, '_voxnormed');
                            auc_voxnormed = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(conditionName).(feat_voxnormed).AUC});
                            height_voxnormed = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(conditionName).(feat_voxnormed).height});
                            time_to_peak_voxnormed = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(conditionName).(feat_voxnormed).time_to_peak});
                            half_height_voxnormed = cell2mat({avgHRF.fit{fitIndex}{regionIndex}.(conditionName).(feat_voxnormed).half_height});
                            
                            avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phase(p).auc = unique(auc(idx));
                            avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phase(p).auc_voxnormed = unique(auc_voxnormed(idx));
    
                            avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phase(p).height = height(idx);
                            avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phase(p).height_voxnormed = height_voxnormed(idx);
                            avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phase(p).time_to_peak = time_to_peak(idx);
                            avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phase(p).time_to_peak_voxnormed = time_to_peak_voxnormed(idx);
                            avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phase(p).half_height = half_height(idx);
                            avgHRF.fit{fitIndex}{regionIndex}.(conditionName).phase(p).half_height_voxnormed = half_height_voxnormed(idx);
                        end
                    end
                end
            end
        end
        % Calculate the within-subject SE for each time point NOT DONE
        % within_subject_SE= squeeze(std(avgVector_across_conditions, 0, 1)) / sqrt(size(models, 1));

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



% Plotting Demonstration:
% 
% % Your time vector, average vector, and standard deviation vector
% time = 1:width(avgVector);  % Replace with your actual time vector
% c=3
% 
% % Create a new figure
% figure;
% 
% % Plot the average vector
% plot(time, avgVector_across_conditions(c,:), 'b', 'LineWidth', 2); 
% hold on;
% 
% % Calculate the upper and lower bounds of the standard deviation
% % upperBound = avgVector + stdVector;
% % lowerBound = avgVector - stdVector;
% 
% % SE
% % upperBound = avgVector_across_conditions(c,:) + stdVector/sqrt(10);
% % lowerBound = avgVector_across_conditions(c,:) - stdVector/sqrt(10);
% 
% 
% % Create a shaded area for the standard deviation
% x = [time, fliplr(time)];
% y = [upperBound, fliplr(lowerBound)];
% fill(x, y, 'b', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
% 
% % Labels and title for clarity
% xlabel('Time');
% ylabel('Measurement');
% title('Average Measurement with Standard Deviation');
% legend('Average', '1 SD', 'Location', 'Best');
% 
% % Display the plot
% hold off;
% grid on;
% 
% 
% 
% 
% 
% % Computing the within-subject SE
% % Assuming `data` is a 3D matrix: subjects x time points x conditions.
% 
% num_subjects=10
% num_timepoints=64
% num_conditions=3
% 
% 
% data=models3d;
% % Assuming `data` is a 3D matrix: subjects x time points x conditions.
% 
% [num_subjects, num_timepoints, num_conditions] = size(data);
% 
% % Calculating the mean across conditions for each subject and timepoint
% subject_means = mean(data, 3); 
% 
% % Calculating the standard deviation of the subject means across timepoints
% subject_std = std(subject_means, 0, 2);
% 
% % Calculating the average standard deviation across subjects
% avg_std = mean(subject_std);
% 
% % Calculating the within-subject standard error
% WSE = avg_std / sqrt(num_conditions);
% 
% % Preallocate arrays for efficiency.
% mean_response = zeros(num_conditions, num_timepoints);
% 
% % Loop through each condition, time point and calculate the mean response.
% for c = 1:num_conditions
%     for t = 1:num_timepoints
%         mean_response(c, t) = mean(data(:, t, c));
%     end
% 
%     % Optionally, plot the mean response with shaded WSE for each condition.
%     figure;
%     plot(mean_response(c, :), 'LineWidth', 2); % Plot the mean response.
%     hold on;
% 
%     % Add shaded error (WSE).
%     fill([1:num_timepoints, fliplr(1:num_timepoints)], ...
%          [mean_response(c, :) + WSE, fliplr(mean_response(c, :) - WSE)], ...
%          [0.9, 0.9, 0.9], 'EdgeColor', 'none');
% 
%     hold off;
%     xlabel('Time Point');
%     ylabel('Response');
%     title(['Mean Response Over Time with WSE for Condition ', num2str(c)]);
% end
% 
% 
% 
% barplot_get_within_ste(models3d(:,:,1))
% barplot_get_within_ste(models3d(:,:,2))
% barplot_get_within_ste(models3d(:,:,3))

