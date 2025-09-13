function [HRF, tc] = extractHRF(HRF_OBJ, CondNames, varargin)
    % Passed in HRF_OBJ is a cell array of fmri_data for each condition.
    % CondNames should be a cell array of charstr for each condition

    isSigFound = false;
    isAtlasFound = false;
    isRegsFound = false;
    r = [];
    s = [];

    if isa(HRF_OBJ, 'fmri_data')
        HRF_OBJ={HRF_OBJ};
    end
    

    % Initialize the parallel pool if it's not already running
    % if isempty(gcp('nocreate'))
    %     parpool;
    % end

    for k = 1:length(varargin)
        if strcmpi(varargin{k}, 'atlas')
            if isa(varargin{k+1}, 'atlas')
                at=varargin{k+1};
            else
                error('Passed in atlas not an atlas.');
            end

            isAtlasFound = true;
        end

        if strcmpi(varargin{k}, 'regions')
            if isAtlasFound
                if ischar(varargin{k+1})
                    r={varargin{k+1}};
                    isRegsFound = true;
                elseif iscell(varargin{k+1})
                    r=varargin{k+1};
                    isRegsFound = true;
                end
            else
                error('Cannot extract HRF of regions without atlas')
            end
        end

        if isAtlasFound && ~isRegsFound
            % Extract all atlas regions.
            r=at.labels;
        end

        if strcmpi(varargin{k}, 'sig') % Expected input: 'nps', 'siips', or 'all'
            if ischar(varargin{k+1})
                s={varargin{k+1}};
                isSigFound = true;
            elseif iscell(varargin{k+1})
                s=varargin{k+1};
                isSigFound = true;
            else
                error(['Input argument for sig is unknown.']);
            end

        end
    end

    % ROIs will consist of regions and signatures
    rois = [r, s];

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
        if isSigFound & strcmpi(rois{r}, 'all')
            [sig_dp, ~]=apply_all_signatures(HRF_OBJ, 'similarity_metric', 'dot_product', 'conditionnames', CondNames, 'image_set', 'all');
            [sig_cs, ~]=apply_all_signatures(HRF_OBJ, 'similarity_metric', 'cosine_similarity', 'conditionnames', CondNames, 'image_set', 'all');
            [sig_r,  ~]=apply_all_signatures(HRF_OBJ, 'similarity_metric', 'correlation', 'conditionnames', CondNames, 'image_set', 'all');
        end


        for c=1:numel(CondNames)

            % Check to see if its a signature name first
            if isSigFound
                
                switch rois{r}
                    case 'all'
                        for sig = 1:numel(sig_dp.signaturenames)
                            try
                                HRF{r}.([CondNames{c},'_',sig_dp.signaturenames{sig},'_dp']).model=sig_dp.(sig_dp.signaturenames{sig}).(CondNames{c});
                                HRF{r}.([CondNames{c},'_',sig_cs.signaturenames{sig},'_cs']).model=sig_cs.(sig_cs.signaturenames{sig}).(CondNames{c});
                                HRF{r}.([CondNames{c},'_',sig_r.signaturenames{sig},'_r']).model=sig_r.(sig_r.signaturenames{sig}).(CondNames{c});
    
                                [HRF{r}.([CondNames{c},'_',sig_dp.signaturenames{sig},'_dp']).peaks, HRF{r}.([CondNames{c},'_',sig_dp.signaturenames{sig},'_dp']).troughs]=detectPeaksTroughs(sig_dp.(sig_dp.signaturenames{sig}).(CondNames{c}), false);
                                [HRF{r}.([CondNames{c},'_',sig_cs.signaturenames{sig},'_cs']).peaks, HRF{r}.([CondNames{c},'_',sig_cs.signaturenames{sig},'_cs']).troughs]=detectPeaksTroughs(sig_cs.(sig_cs.signaturenames{sig}).(CondNames{c}), false);
                                [HRF{r}.([CondNames{c},'_',sig_r.signaturenames{sig},'_r']).peaks, HRF{r}.([CondNames{c},'_',sig_r.signaturenames{sig},'_r']).troughs]=detectPeaksTroughs(sig_r.(sig_r.signaturenames{sig}).(CondNames{c}), false);
    
                            catch
                                disp(c);
                                disp(rois{r});
                            end
                        end

                    case 'nps'
                        % tc{r}=apply_nps(HRF_OBJ); % This does all conditions at once!
                        try
                            % tc{r}{c}=apply_nps(HRF_OBJ{c});
                            HRF{r}.([CondNames{c}, '_NPS_dp']).model=cell2mat(apply_nps(HRF_OBJ{c}));
                            HRF{r}.([CondNames{c}, '_NPS_cs']).model=cell2mat(apply_nps(HRF_OBJ{c}, 'cosine_similarity'));
                            HRF{r}.([CondNames{c}, '_NPS_r']).model=cell2mat(apply_nps(HRF_OBJ{c}, 'correlation'));
                            [HRF{r}.([CondNames{c}, '_NPS_dp']).peaks, HRF{r}.([CondNames{c}, '_NPS_dp']).troughs]=detectPeaksTroughs(HRF{r}.([CondNames{c}, '_NPS_dp']).model, false);
                            [HRF{r}.([CondNames{c}, '_NPS_cs']).peaks, HRF{r}.([CondNames{c}, '_NPS_cs']).troughs]=detectPeaksTroughs( HRF{r}.([CondNames{c}, '_NPS_cs']).model, false);
                            [HRF{r}.([CondNames{c}, '_NPS_r']).peaks, HRF{r}.([CondNames{c}, '_NPS_r']).troughs]=detectPeaksTroughs(HRF{r}.([CondNames{c}, '_NPS_r']).model, false);
                        catch
                            disp(c);
                            disp(rois{r});
                        end

                    case 'siips'
                        try
                            HRF{r}.([CondNames{c}, '_SIIPS_dp']).model=cell2mat(apply_siips(HRF_OBJ{c}));
                            HRF{r}.([CondNames{c}, '_SIIPS_cs']).model=cell2mat(apply_siips(HRF_OBJ{c}, 'cosine_similarity'));
                            HRF{r}.([CondNames{c}, '_SIIPS_r']).model=cell2mat(apply_siips(HRF_OBJ{c}, 'correlation'));
                            [HRF{r}.([CondNames{c}, '_SIIPS_dp']).peaks, HRF{r}.([CondNames{c}, '_SIIPS_dp']).troughs]=detectPeaksTroughs(HRF{r}.([CondNames{c}, '_SIIPS_dp']).model, false);
                            [HRF{r}.([CondNames{c}, '_SIIPS_cs']).peaks, HRF{r}.([CondNames{c}, '_SIIPS_cs']).troughs]=detectPeaksTroughs( HRF{r}.([CondNames{c}, '_SIIPS_cs']).model, false);
                            [HRF{r}.([CondNames{c}, '_SIIPS_r']).peaks, HRF{r}.([CondNames{c}, '_SIIPS_r']).troughs]=detectPeaksTroughs(HRF{r}.([CondNames{c}, '_SIIPS_r']).model, false);
                        catch
                            disp(c);
                            disp(rois{r});
                        end
                end
            end


            if isAtlasFound & ~ismember(rois{r}, {'all', 'nps', 'siips'})


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
            
    end

    % HRF_struct=struct;
    % if isAtlasFound
    %     HRF_struct.atlas=at;
    % end
    % HRF_struct.region=rois;
    % HRF_struct.CondNames=CondNames;
    % HRF_struct.fit=HRF;
    % HRF=HRF_struct;
    


    % temp_HRF_fit = HRF_local;
    % Save the results for this ROI
    % display([num2str(t), ' Done!'])
    % temp_HRF_fit{t} = HRF_local;
    % tc{t} = tc;


    % Transfer the results from the temporary cell array to the HRF structure
    % HRF.fit = temp_HRF_fit;
    % HRF.params=HRF_PARAMS;

    % delete(gcp('nocreate'));



end