function [tc, HRF]=EstimateHRF_inAtlas(fmri_d, PREPROC_PARAMS, HRF_PARAMS, at, rois, outfile)
    % EstimateHRF_inAtlas takes a raw 4D fmri_data object, preprocesses it, and
    % then outputs an estimated HRF time series for each condition of interest.
    % PREPROC struct needs to have TR, hpf, and Condition information
    % HRF struct is a structure of HRF fitting parameters. It needs to have T, FWHM, alpha, and type, otherwise default
    % values will be supplied.
    %
    % Michael Sun, Ph.D.
    % - Takes 4D fmri_data() object
    % - PREPROC_PARAMS: struct object that contains the fields: TR, R, hpf, and smooth
    % - HRF_PARAMS: struct object that contains the fields: Condition, CondNames, TR, T, FWHM, alpha, and types. 
    % - at: atlas object
    % - rois: cell-array of labels that match labels in at.labels
    % - outfile: Desired filepath without extension. Extensions will be appended.
    %
    % *Usage:
    % ::
    %    [tc, HRF] = EstimateHRF_inAtlas(image_obj, PREPROC_PARAMS, HRF_PARAMS, at, rois, outfile})
    %
    % TODO: Pass in estHRF directory to reuse nifti files.

    % Step 0. Check if a directory exists that has estHRF .nii outputs
    % already.


    % Step 1. Preprocess 
    [preproc_dat]=canlab_connectivity_preproc(fmri_d, PREPROC_PARAMS.R, 'hpf', PREPROC_PARAMS.hpf, PREPROC_PARAMS.TR, 'average_over', 'no_plots');
    
    % Step 2. Smooth
    preproc_dat=preprocess(preproc_dat, 'smooth', PREPROC_PARAMS.smooth);
    
    % Step 3. fitHRF to each ROI's worth of data
    [tc, HRF]=roiTS_fitHRF(preproc_dat, HRF_PARAMS, rois, at, outfile);
    HRF.preproc_params=PREPROC_PARAMS;
    
    % Step 4. Save your Results
    try
        % saveas(gcf, [outfile, '.png'])
        save([outfile, '.mat'], 'tc', 'HRF', '-v7.3');
    catch
        warning([outfile, ' files could not be saved.']);
    end
end


%% HELPER FUNCTIONS

function [tc, HRF]=roiTS_fitHRF(preproc_dat, HRF_PARAMS, rois, at, outfile)

    % Make directories for files if needed
    if ~isempty(fileparts(outfile))
        if ~exist(fileparts(outfile), 'dir')
            mkdir(fileparts(outfile));
        end
    end

    HRF.atlas=at.atlas_name;
    HRF.region=rois;
    HRF.types=HRF_PARAMS.types;
    HRF.name=preproc_dat.image_names;
    
    % Initialize the parallel pool if it's not already running
    if isempty(gcp('nocreate'))
        parpool;
    end

    % Initialize 'tc' and 'temp_HRF_fit' cell arrays
    tc = cell(1, numel(HRF_PARAMS.types));
    temp_HRF_fit = cell(1, numel(HRF_PARAMS.types));

    % Write out the images for later post-analyses
    % [~, fname, ~]=fileparts(preproc_dat.image_names);
    % fname=outfile

    parfor t=1:numel(HRF_PARAMS.types)
    % for t=1:numel(HRF_PARAMS.types)  % FOR TROUBLESHOOTING  
        warning('off', 'all');
        if strcmp(HRF_PARAMS.types{t}, 'IL')

            [~, ~, PARAM_OBJ, HRF_OBJ] = hrf_fit(preproc_dat, HRF_PARAMS.TR, HRF_PARAMS.Condition, HRF_PARAMS.T, HRF_PARAMS.types{t}, 0);

        elseif strcmp(HRF_PARAMS.types{t}, 'FIR')

            [~, ~, PARAM_OBJ, HRF_OBJ] = hrf_fit(preproc_dat, HRF_PARAMS.TR, HRF_PARAMS.Condition, HRF_PARAMS.T, HRF_PARAMS.types{t}, 1);


        elseif strcmp(HRF_PARAMS.types{t}, 'CHRF')

            [~, ~, PARAM_OBJ, HRF_OBJ] = hrf_fit(preproc_dat, HRF_PARAMS.TR, HRF_PARAMS.Condition, HRF_PARAMS.T, HRF_PARAMS.types{t}, 2);
        else
            error('No valid fit-type. Choose IL, FIR, or CHRF')
        end
        
        for c=1:numel(HRF_PARAMS.Condition)

            HRF_OBJ{c}.fullpath=sprintf([outfile, '_type-', HRF_PARAMS.types{t}, '_condition-', HRF_PARAMS.CondNames{c}, '_fit.nii']);
            PARAM_OBJ{c}.fullpath=sprintf([outfile, '_type-', HRF_PARAMS.types{t}, '_condition-', HRF_PARAMS.CondNames{c}, '_params.nii']);
            try
                write(HRF_OBJ{c}, 'overwrite');
                write(PARAM_OBJ{c}, 'overwrite');
            catch
                warning('Not able to write one or more files.');
            end

        end

        HRF_local = cell(1, numel(HRF_PARAMS.types));
        tc_local = cell(1, numel(rois));


        % Consider doing apply_parcellation instead of mean(apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat);
        % [parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg] = apply_parcellation(dat,at);
        % nps=load_image_set('npsplus');
        % nps = get_wh_image(nps,1);
        % [parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg] = apply_parcellation(dat,at, 'pattern_expression', nps);
        % r=region(at,'unique_mask_values');
        % wh_parcels=~all(isnan(parcel_means))
        
        for r=1:numel(rois)
            tic
            for c=1:numel(HRF_PARAMS.CondNames)

                try
                    tc_local{r}{c}=mean(apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat);
            
                    HRF_local{r}.(HRF_PARAMS.CondNames{c}).model=tc_local{r}{c};
                    [HRF_local{r}.(HRF_PARAMS.CondNames{c}).peaks, HRF_local{r}.(HRF_PARAMS.CondNames{c}).troughs]=detectPeaksTroughs(tc_local{r}{c}', false);

                    [~, regionVoxNum, ~, ~]=at.select_atlas_subset(rois(r), 'exact').get_region_volumes;
                    HRF_local{r}.(HRF_PARAMS.CondNames{c}).model_voxnormed=tc_local{r}{c}/regionVoxNum;
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



