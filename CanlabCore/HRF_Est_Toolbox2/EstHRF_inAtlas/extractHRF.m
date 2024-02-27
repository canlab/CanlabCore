function [HRF, tc] = extractHRF(HRF_OBJ, CondNames, at, rois)
    % Passed in HRF_OBJ is a cell array of fmri_data for each condition.
    % CondNames should be a cell array of charstr for each condition


    % Initialize the parallel pool if it's not already running
    if isempty(gcp('nocreate'))
        parpool;
    end

    % Preallocation of tc_local before the parfor loop
    HRF= cell(1, numel(rois));
    tc= cell(1, numel(rois));

    CondNames=matlab.lang.makeValidName(CondNames);

    % Now, handle the fourth dimension which varies with 'd'
    numCondNames = numel(CondNames);

    % HRF_local{d} = cell(1, numel(rois));
    % tc_local{d} = cell(1, numel(rois));

    % Consider doing apply_parcellation instead of mean(apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat);
    % [parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg] = apply_parcellation(dat,at);
    % nps=load_image_set('npsplus');
    % nps = get_wh_image(nps,1);
    % [parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg] = apply_parcellation(dat,at, 'pattern_expression', nps);
    % r=region(at,'unique_mask_values');
    % wh_parcels=~all(isnan(parcel_means))
    
    for r=1:numel(rois)
        HRF{r} = struct;
        tc{r} = cell(1, numCondNames);
        tic
        for c=1:numel(CondNames)

            try
                tc{r}{c}=mean(apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat);
        
                HRF{r}.(CondNames{c}).model=tc{r}{c};
                [HRF{r}.(CondNames{c}).peaks, HRF{r}.(CondNames{c}).troughs]=detectPeaksTroughs(tc{r}{c}', false);

                [~, regionVoxNum, ~, ~]=at.select_atlas_subset(rois(r), 'exact').get_region_volumes;
                HRF{r}.(CondNames{c}).model_voxnormed=tc{r}{c}/regionVoxNum;
                [HRF{r}.(CondNames{c}).peaks_voxnormed, HRF{r}.(CondNames{c}).troughs_voxnormed]=detectPeaksTroughs(tc{r}{c}'/regionVoxNum, false);
            catch
                % disp(t);
                disp(c);
                disp(rois{r});
                tc{r}{c}
                mean(apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat)
                {apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat}
                HRF_OBJ

            end

            % Number of phases
            start_times = [HRF{r}.(CondNames{c}).peaks.start_time, HRF{r}.(CondNames{c}).troughs.start_time];
            end_times = [HRF{r}.(CondNames{c}).peaks.end_time, HRF{r}.(CondNames{c}).troughs.end_time];
            phases = [start_times(:), end_times(:)];
            unique_phases = unique(phases, 'rows');
            HRF{r}.(CondNames{c}).phases = mat2cell(unique_phases, ones(size(unique_phases, 1), 1), 2);

            % Loop over phases
            for p = 1:numel(HRF{r}.(CondNames{c}).phases)
                features = {'peaks', 'troughs'};
                for f = 1:2
                    feat = features{f};
                    
                    % Count the number of features (peaks or troughs)
                    start_times = cell2mat({HRF{r}.(CondNames{c}).(feat).start_time});
                    current_phase_start = HRF{r}.(CondNames{c}).phases{p}(1);
                    HRF{r}.(CondNames{c}).phase(p).(feat) = sum(start_times == current_phase_start);
                    
                    if HRF{r}.(CondNames{c}).phase(p).(feat) > 0
                        display(['Estimating ' , '_', rois{r}, '_', CondNames{c}, '_', ' Now...!']);
                        display(['Phase ' , num2str(p), 'Feature ', feat]);
                        idx = find(start_times == current_phase_start);
                        auc = cell2mat({HRF{r}.(CondNames{c}).(feat).AUC});
                        height = cell2mat({HRF{r}.(CondNames{c}).(feat).height});
                        time_to_peak = cell2mat({HRF{r}.(CondNames{c}).(feat).time_to_peak});
                        half_height = cell2mat({HRF{r}.(CondNames{c}).(feat).half_height});
                        
                        feat_voxnormed = strcat(feat, '_voxnormed');
                        auc_voxnormed = cell2mat({HRF{r}.(CondNames{c}).(feat_voxnormed).AUC});
                        height_voxnormed = cell2mat({HRF{r}.(CondNames{c}).(feat_voxnormed).height});
                        time_to_peak_voxnormed = cell2mat({HRF{r}.(CondNames{c}).(feat_voxnormed).time_to_peak});
                        half_height_voxnormed = cell2mat({HRF{r}.(CondNames{c}).(feat_voxnormed).half_height});
                        
                        HRF{r}.(CondNames{c}).phase(p).auc = unique(auc(idx));
                        HRF{r}.(CondNames{c}).phase(p).auc_voxnormed = unique(auc_voxnormed(idx));

                        HRF{r}.(CondNames{c}).phase(p).height = height(idx);
                        HRF{r}.(CondNames{c}).phase(p).height_voxnormed = height_voxnormed(idx);
                        HRF{r}.(CondNames{c}).phase(p).time_to_peak = time_to_peak(idx);
                        HRF{r}.(CondNames{c}).phase(p).time_to_peak_voxnormed = time_to_peak_voxnormed(idx);
                        HRF{r}.(CondNames{c}).phase(p).half_height = half_height(idx);
                        HRF{r}.(CondNames{c}).phase(p).half_height_voxnormed = half_height_voxnormed(idx);
                    end
                end
            end
        end
        display(strjoin({' Done in ', num2str(toc), ' seconds with ', rois{r}})); 
        
    end


    % temp_HRF_fit = HRF_local;
    % Save the results for this ROI
    % display([num2str(t), ' Done!'])
    % temp_HRF_fit{t} = HRF_local;
    % tc{t} = tc;


    % Transfer the results from the temporary cell array to the HRF structure
    % HRF.fit = temp_HRF_fit;
    % HRF.params=HRF_PARAMS;

    delete(gcp('nocreate'));



end