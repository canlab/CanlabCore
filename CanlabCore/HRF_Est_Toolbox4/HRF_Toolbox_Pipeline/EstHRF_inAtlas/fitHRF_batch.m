function [tc, HRF]=fitHRF_batch(datadir, HRF_PARAMS, rois, at, outfile)
    % This is not done yet.

    % Make directories for files if needed
    if ~isempty(fileparts(outfile))
        if ~exist(fileparts(outfile), 'dir')
            mkdir(fileparts(outfile));
        end
    end

    HRF.atlas=at.atlas_name;
    HRF.region=rois;
    HRF.types=HRF_PARAMS.types;
    HRF.name=fmri_d.image_names;
    
    % Initialize 'tc' and 'temp_HRF_fit' cell arrays
    tc = cell(1, numel(HRF_PARAMS.types);
    temp_HRF_fit = cell(1, numel(HRF_PARAMS.types));

    % Write out the images for later post-analyses
    % [~, fname, ~]=fileparts(preproc_dat.image_names);
    % fname=outfile

    HRF_local = cell(1, numel(HRF_PARAMS.types));
    tc_local = cell(1, numel(rois));

    
    canlab_list_subjects()

    fmri_data()



    for r=1:numel(rois)
        tic
        for c=1:numel(HRF_PARAMS.CondNames)

            try
                tc_local{r}{c}=mean(apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat);
        
                HRF_local{r}.(HRF_PARAMS.CondNames{c}).model=tc_local{r}{c};
                [HRF_local{r}.(HRF_PARAMS.CondNames{c}).peaks, HRF_local{r}.(HRF_PARAMS.CondNames{c}).troughs]=detectPeaksTroughs(tc_local{r}{c}', false);

                [~, regionVoxNum, ~, ~]=at.select_atlas_subset(rois(r), 'exact').get_region_volumes;
                HRF_local{r}.(HRF_PARAMS.CondNames{c}).model_voxnormed=tc_local{r}/regionVoxNum;
                [HRF_local{r}.(HRF_PARAMS.CondNames{c}).peaks_voxnormed, HRF_local{r}.(HRF_PARAMS.CondNames{c}).troughs_voxnormed]=detectPeaksTroughs(tc_local{r}{c}'/regionVoxNum, false);
            catch
                disp(t);
                disp(c);
                disp(rois{r});
                tc_local{r}{c}
                mean(apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat)
                {apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat}
                HRF_OBJ{c}

            end

            % Number of phases
            start_times = [HRF_local{r}.(HRF_PARAMS.CondNames{c}).peaks.start_time, HRF_local{r}.(HRF_PARAMS.CondNames{c}).troughs.start_time];
            end_times = [HRF_local{r}.(HRF_PARAMS.CondNames{c}).peaks.end_time, HRF_local{r}.(HRF_PARAMS.CondNames{c}).troughs.end_time];
            phases = [start_times(:), end_times(:)];
            unique_phases = unique(phases, 'rows');
            HRF_local{r}.(HRF_PARAMS.CondNames{c}).phases = mat2cell(unique_phases, ones(size(unique_phases, 1), 1), 2);

            % Loop over phases
            for p = 1:numel(HRF_local{r}.(HRF_PARAMS.CondNames{c}).phases)
                features = {'peaks', 'troughs'};
                for f = 1:2
                    feat = features{f};
                    
                    % Count the number of features (peaks or troughs)
                    start_times = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat).start_time});
                    current_phase_start = HRF_local{r}.(HRF_PARAMS.CondNames{c}).phases{p}(1);
                    HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).(feat) = sum(start_times == current_phase_start);
                    
                    if HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).(feat) > 0
                        display(['Estimating ' , num2str(t), '_', rois{r}, '_', HRF_PARAMS.CondNames{c}, '_', ' Now...!'])
                        display(['Phase ' , num2str(p), 'Feature ', feat])
                        idx = find(start_times == current_phase_start)
                        auc = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat).AUC})
                        height = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat).height})
                        time_to_peak = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat).time_to_peak})
                        half_height = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat).half_height})
                        
                        feat_voxnormed = strcat(feat, '_voxnormed')
                        auc_voxnormed = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat_voxnormed).AUC})
                        height_voxnormed = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat_voxnormed).height})
                        time_to_peak_voxnormed = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat_voxnormed).time_to_peak})
                        half_height_voxnormed = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat_voxnormed).half_height})
                        
                        HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).auc = unique(auc(idx))
                        HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).auc_voxnormed = unique(auc_voxnormed(idx))

                        HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).height = height(idx)
                        HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).height_voxnormed = height_voxnormed(idx)
                        HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).time_to_peak = time_to_peak(idx)
                        HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).time_to_peak_voxnormed = time_to_peak_voxnormed(idx)
                        HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).half_height = half_height(idx)
                        HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).half_height_voxnormed = half_height_voxnormed(idx)
                    end
                end
            end
        end
        display([num2str(t), ' Done in ', toc, ' with ' rois(r)]); 
        
    end
    % Save the results for this ROI
    display([num2str(t), ' Done!'])
    temp_HRF_fit{t} = HRF_local;
    tc{t} = tc_local;
end

% Transfer the results from the temporary cell array to the HRF structure
HRF.fit = temp_HRF_fit;
HRF.params=HRF_PARAMS;
end