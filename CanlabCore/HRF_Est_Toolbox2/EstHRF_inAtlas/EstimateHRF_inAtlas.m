function [tc, HRF, HRF_OBJ, PARAM_OBJ]=EstimateHRF_inAtlas(fmri_d, PREPROC_PARAMS, HRF_PARAMS, at, rois, outfile)
    % EstimateHRF_inAtlas takes a raw 4D fmri_data object, preprocesses it, and
    % then outputs an estimated HRF time series for each condition of interest.
    % PREPROC struct needs to have TR, hpf, and Condition information
    % HRF struct is a structure of HRF fitting parameters. It needs to have T, FWHM, alpha, and type, otherwise default
    % values will be supplied.
    %
    % Michael Sun, Ph.D.
    % - Takes 4D fmri_data() object or cell-array of fmri_data() 
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

    % Step 0. Load in an SPM.mat file if available to pass in metadata and
    % such


    % Check if its an SPM struct or filepath:
    if ischar(fmri_d)
        if contains(fmri_d, 'SPM.mat')
            load(fmri_d);
            if exist('SPM', 'var')
                isSPM=true;
            end
        
        else
            fmri_d=fmri_data(fmri_d);
        end
    end
    
    if isstruct(fmri_d)
        SPM=fmri_d;
        isSPM=true;
    end



    if isSPM
        % if isstruct(varargin{1})
        %     SPM=varargin{1};
        % else
        %     % If pointing to an SPM filepath
        %     load(varargin{1})
        % end
        % % Don't do that, spmify instead into gKWY (grand-mean-scaled, filtered (K), and whitened (W))
        % gKWY=spmify(fmri_d, SPM);
        % 
        % for d=1:numel(gKWY)
        %     preproc_dat{d}=fmri_d{d};
        %     preproc_dat{d}.dat=gKWY{d};
        % end
        [tc, HRF, HRF_OBJ, PARAM_OBJ]=roiTS_fitHRF_SPM(SPM, HRF_PARAMS, rois, at, outfile);

    else
        % Step 1. Preprocess
        for d=1:numel(fmri_d)
            [preproc_dat{d}]=canlab_connectivity_preproc(fmri_d, PREPROC_PARAMS.R{d}, 'hpf', PREPROC_PARAMS.hpf, PREPROC_PARAMS.TR, 'average_over', 'no_plots');
            % Step 2. Smooth
            preproc_dat{d}=preprocess(preproc_dat{d}, 'smooth', PREPROC_PARAMS.smooth);
        end
    
        % Step 3. fitHRF to each ROI's worth of data
        [tc, HRF]=roiTS_fitHRF(preproc_dat, HRF_PARAMS, rois, at, outfile);


    end


    if numel(HRF)>1
        for i=1:numel(HRF)
            HRF{i}.preproc_params=PREPROC_PARAMS;        
        end    
    else
        HRF.preproc_params=PREPROC_PARAMS;
    end
    
    
    % Step 4. Save your Results
    try
        % saveas(gcf, [outfile, '.png'])
        save([outfile, '.mat'], 'tc', 'HRF', '-v7.3');
    catch
        warning([outfile, ' files could not be saved.']);
    end
end


%% HELPER FUNCTIONS
% function [HRF_OBJ, HRF]=gen_HRFimg(preproc_dat, HRF_PARAMS, rois, at, outfile)
% 
%     % This is finished but untested. Intention is to generate a directory from
%     % HRF images in BIDS format.
% 
%     % Make directories for files if needed
%     if ~isempty(fileparts(outfile))
%         if ~exist(fileparts(outfile), 'dir')
%             mkdir(fileparts(outfile));
%         end
%     end
% 
%     % Make sure CondNames are valid before continuing:
%     HRF_PARAMS.CondNames=matlab.lang.makeValidName(HRF_PARAMS.CondNames);
% 
%     HRF.atlas=at;
%     HRF.region=rois;
%     HRF.types=HRF_PARAMS.types;
%     HRF.name=preproc_dat.image_names;
% 
%     % Initialize the parallel pool if it's not already running
%     if isempty(gcp('nocreate'))
%         parpool;
%     end
% 
%     % Write out the images for later post-analyses
%     % [~, fname, ~]=fileparts(preproc_dat.image_names);
%     % fname=outfile
% 
%     parfor t=1:numel(HRF_PARAMS.types)
%     % for t=1:numel(HRF_PARAMS.types)  % FOR TROUBLESHOOTING  
%         warning('off', 'all');
%         switch HRF_PARAMS.types{t}
%             case 'IL'
%                 [~, ~, PARAM_OBJ, HRF_OBJ] = hrf_fit(preproc_dat, HRF_PARAMS.TR, HRF_PARAMS.Condition, HRF_PARAMS.T, HRF_PARAMS.types{t}, 0);
% 
%             case 'FIR'
%                 [~, ~, PARAM_OBJ, HRF_OBJ] = hrf_fit(preproc_dat, HRF_PARAMS.TR, HRF_PARAMS.Condition, HRF_PARAMS.T, HRF_PARAMS.types{t}, 1);
% 
%             case 'CHRF'
%                 [~, ~, PARAM_OBJ, HRF_OBJ] = hrf_fit(preproc_dat, HRF_PARAMS.TR, HRF_PARAMS.Condition, HRF_PARAMS.T, HRF_PARAMS.types{t}, 2);
%             otherwise
%                 error('No valid fit-type. Choose IL, FIR, or CHRF')
%         end
% 
%         for c=1:numel(HRF_PARAMS.Condition)
% 
%             HRF_OBJ{c}.fullpath=sprintf([outfile, '_type-', HRF_PARAMS.types{t}, '_condition-', HRF_PARAMS.CondNames{c}, '_fit.nii']);
%             PARAM_OBJ{c}.fullpath=sprintf([outfile, '_type-', HRF_PARAMS.types{t}, '_condition-', HRF_PARAMS.CondNames{c}, '_params.nii']);
%             try
%                 write(HRF_OBJ{c}, 'overwrite');
%                 write(PARAM_OBJ{c}, 'overwrite');
%             catch
%                 warning('Not able to write one or more files.');
%             end
% 
%         end
% 
%     end
% end


% function [tc, HRF]=gen_HRFstruct_from_dir(preproc_dat, HRF_PARAMS, rois, at, outfile, estHRF_dir, HRF_OBJ)
% 
%     % This is unfinished. Intention is to take a generated directory from
%     % gen_HRFimg and generate an HRF structure from it.
% 
%     % Initialize the parallel pool if it's not already running
%     if isempty(gcp('nocreate'))
%         parpool;
%     end
% 
%     % Initialize 'tc' and 'temp_HRF_fit' cell arrays
%     tc = cell(1, numel(HRF_PARAMS.types));
%     temp_HRF_fit = cell(1, numel(HRF_PARAMS.types));
% 
%     parfor t=1:numel(HRF_PARAMS.types)
% 
%         HRF_local = cell(1, numel(rois));
%         tc_local = cell(1, numel(rois));
% 
% 
%         % Consider doing apply_parcellation instead of mean(apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat);
%         % [parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg] = apply_parcellation(dat,at);
%         % nps=load_image_set('npsplus');
%         % nps = get_wh_image(nps,1);
%         % [parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg] = apply_parcellation(dat,at, 'pattern_expression', nps);
%         % r=region(at,'unique_mask_values');
%         % wh_parcels=~all(isnan(parcel_means))
% 
%         for r=1:numel(rois)
%             tic
%             for c=1:numel(HRF_PARAMS.CondNames)
% 
%                 try
%                     tc_local{r}{c}=mean(apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat);
% 
%                     HRF_local{r}.(HRF_PARAMS.CondNames{c}).model=tc_local{r}{c};
%                     [HRF_local{r}.(HRF_PARAMS.CondNames{c}).peaks, HRF_local{r}.(HRF_PARAMS.CondNames{c}).troughs]=detectPeaksTroughs(tc_local{r}{c}', false);
% 
%                     [~, regionVoxNum, ~, ~]=at.select_atlas_subset(rois(r), 'exact').get_region_volumes;
%                     HRF_local{r}.(HRF_PARAMS.CondNames{c}).model_voxnormed=tc_local{r}{c}/regionVoxNum;
%                     [HRF_local{r}.(HRF_PARAMS.CondNames{c}).peaks_voxnormed, HRF_local{r}.(HRF_PARAMS.CondNames{c}).troughs_voxnormed]=detectPeaksTroughs(tc_local{r}{c}'/regionVoxNum, false);
%                 catch
%                     disp(t);
%                     disp(c);
%                     disp(rois{r});
%                     tc_local{r}{c}
%                     mean(apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat)
%                     {apply_mask(HRF_OBJ{c}, at.select_atlas_subset(rois(r), 'exact')).dat}
%                     HRF_OBJ{c}
% 
%                 end
% 
%                 % Number of phases
%                 start_times = [HRF_local{r}.(HRF_PARAMS.CondNames{c}).peaks.start_time, HRF_local{r}.(HRF_PARAMS.CondNames{c}).troughs.start_time];
%                 end_times = [HRF_local{r}.(HRF_PARAMS.CondNames{c}).peaks.end_time, HRF_local{r}.(HRF_PARAMS.CondNames{c}).troughs.end_time];
%                 phases = [start_times(:), end_times(:)];
%                 unique_phases = unique(phases, 'rows');
%                 HRF_local{r}.(HRF_PARAMS.CondNames{c}).phases = mat2cell(unique_phases, ones(size(unique_phases, 1), 1), 2);
% 
%                 % Loop over phases
%                 for p = 1:numel(HRF_local{r}.(HRF_PARAMS.CondNames{c}).phases)
%                     features = {'peaks', 'troughs'};
%                     for f = 1:2
%                         feat = features{f};
% 
%                         % Count the number of features (peaks or troughs)
%                         start_times = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat).start_time});
%                         current_phase_start = HRF_local{r}.(HRF_PARAMS.CondNames{c}).phases{p}(1);
%                         HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).(feat) = sum(start_times == current_phase_start);
% 
%                         if HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).(feat) > 0
%                             display(['Estimating ' , num2str(t), '_', rois{r}, '_', HRF_PARAMS.CondNames{c}, '_', ' Now...!']);
%                             display(['Phase ' , num2str(p), 'Feature ', feat]);
%                             idx = find(start_times == current_phase_start);
%                             auc = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat).AUC});
%                             height = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat).height});
%                             time_to_peak = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat).time_to_peak});
%                             half_height = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat).half_height});
% 
%                             feat_voxnormed = strcat(feat, '_voxnormed');
%                             auc_voxnormed = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat_voxnormed).AUC});
%                             height_voxnormed = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat_voxnormed).height});
%                             time_to_peak_voxnormed = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat_voxnormed).time_to_peak});
%                             half_height_voxnormed = cell2mat({HRF_local{r}.(HRF_PARAMS.CondNames{c}).(feat_voxnormed).half_height});
% 
%                             HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).auc = unique(auc(idx));
%                             HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).auc_voxnormed = unique(auc_voxnormed(idx));
% 
%                             HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).height = height(idx);
%                             HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).height_voxnormed = height_voxnormed(idx);
%                             HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).time_to_peak = time_to_peak(idx);
%                             HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).time_to_peak_voxnormed = time_to_peak_voxnormed(idx);
%                             HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).half_height = half_height(idx);
%                             HRF_local{r}.(HRF_PARAMS.CondNames{c}).phase(p).half_height_voxnormed = half_height_voxnormed(idx);
%                         end
%                     end
%                 end
%             end
%             display(strjoin({num2str(t), ' Done in ', num2str(toc), ' seconds with ', rois{r}})); 
% 
%         end
%         % Save the results for this ROI
%         display([num2str(t), ' Done!'])
%         temp_HRF_fit{t} = HRF_local;
%         tc{t} = tc_local;
%     end
% 
%     % Transfer the results from the temporary cell array to the HRF structure
%     HRF.fit = temp_HRF_fit;
%     HRF.params=HRF_PARAMS;
% 
%     delete(gcp('nocreate'));
% 
% 
% end


function [tc, HRF]=roiTS_fitHRF(preproc_dat, HRF_PARAMS, rois, at, outfile, varargin)

    % Make directories for files if needed
    if ~isempty(fileparts(outfile))
        disp(fileparts(outfile));
        if ~exist(fileparts(outfile), 'dir')
            mkdir(fileparts(outfile));
        end
    end

    % Make sure CondNames are valid before continuing:
    HRF.atlas=at;
    HRF.region=rois;
    HRF.types=HRF_PARAMS.types;

    if iscell(preproc_dat)
        for c=1:numel(preproc_dat)
            HRF.name{c}=preproc_dat{c}.image_names;
    
            % Initialize 'tc' and 'temp_HRF_fit' cell arrays
            tc{c} = cell(1, numel(HRF_PARAMS.types));
            temp_HRF_fit{c} = cell(1, numel(HRF_PARAMS.types));
        end
    else
        HRF.name=preproc_dat.image_names;

        % Initialize 'tc' and 'temp_HRF_fit' cell arrays
        tc = cell(1, numel(HRF_PARAMS.types));
        temp_HRF_fit = cell(1, numel(HRF_PARAMS.types));
    end

    % Carve up SPM's design matrix for each image
    % if ~isempty(varargin)
    %     if ischar(varargin{1}) || isstring(varargin{1})
    %         load(varargin{1});
    %     elseif isstruct(varargin{1})
    %         SPM=varargin{1};
    %     end

        
    %     DX=cell(1,numel(SPM.nscan));
    %     if numel(preproc_dat) == numel(SPM.nscan)
    %         for d=1:numel(preproc_dat)
    %             % Check if the SPM file accounts for the same number of scans
    % 
    %             % HRF_PARAMS.CondNames
    %             % Use regexp to search for the pattern
    %             % hot_matches = ~cellfun(@isempty, regexp(SPM.xX.name, ['Sn\(', num2str(i), '\).*_heat_start']));
    %             % warm_matches = ~cellfun(@isempty, regexp(SPM.xX.name, ['Sn\(', num2str(i), '\).*_warm_start']));
    %             % imgcue_matches = ~cellfun(@isempty, regexp(SPM.xX.name, ['Sn\(', num2str(i), '\).*_imagine_cue']));
    %             % imag_matches = ~cellfun(@isempty, regexp(SPM.xX.name, ['Sn\(', num2str(i), '\).*_imagine_start']));
    % 
    %             R_matches = ~cellfun(@isempty, regexp(SPM.xX.name, ['Sn\(', num2str(d), '\).*R.*']));
    %             constant_matches = ~cellfun(@isempty, regexp(SPM.xX.name, ['Sn\(', num2str(d), '\).*constant']));
    %             sess_cols = ~cellfun(@isempty, regexp(SPM.xX.name, ['Sn\(', num2str(d), '\).*']));
    % 
    % 
    %             % Find indices of matches
    %             indices = [find(sess_cols)];
    %             % indices = [find(hot_matches) find(warm_matches) find(imgcue_matches) find(imag_matches)];
    %             % indices = [find(R_matches) find(constant_matches)];
    %             task_regressors{d} = [find(sess_cols & (~R_matches & ~constant_matches))];
    % 
    %             % Each design matrix needs task regressors, covariates, and
    %             % intercept
    %             DX{d}=SPM.xX.xKXs.X(SPM.Sess(d).row, indices);
    % 
    %         end
    %         customDX=1;
    %     end
    % 
    % 
    % end

    HRF_PARAMS.CondNames=matlab.lang.makeValidName(HRF_PARAMS.CondNames);



    % Initialize the parallel pool if it's not already running
    if isempty(gcp('nocreate'))
        parpool;
    end


    % Write out the images for later post-analyses
    % [~, fname, ~]=fileparts(preproc_dat.image_names);
    % fname=outfile

    HRF_OBJ=cell(1,numel(preproc_dat));
    PARAM_OBJ=cell(1,numel(preproc_dat));

    parfor d=1:numel(preproc_dat)
    % for d=1:numel(preproc_dat)

        if iscell(preproc_dat)
            data=preproc_dat{d};
        end

        for t=1:numel(HRF_PARAMS.types)


            warning('off', 'all');
            switch HRF_PARAMS.types{t}
                case 'IL'
                    [~, ~, PARAM_OBJ{d}, HRF_OBJ{d}] = hrf_fit(data, HRF_PARAMS.TR, HRF_PARAMS.Condition{d}, HRF_PARAMS.T, HRF_PARAMS.types{t}, 0);
                case 'FIR'
                    [~, ~, PARAM_OBJ{d}, HRF_OBJ{d}] = hrf_fit(data, HRF_PARAMS.TR, HRF_PARAMS.Condition{d}, HRF_PARAMS.T, HRF_PARAMS.types{t}, 0);
                case 'sFIR'
                    [~, ~, PARAM_OBJ{d}, HRF_OBJ{d}] = hrf_fit(data, HRF_PARAMS.TR, HRF_PARAMS.Condition{d}, HRF_PARAMS.T, HRF_PARAMS.types{t}, 1);
                case 'CHRF0'
                    [~, ~, PARAM_OBJ{d}, HRF_OBJ{d}] = hrf_fit(data, HRF_PARAMS.TR, HRF_PARAMS.Condition{d}, HRF_PARAMS.T, HRF_PARAMS.types{t}, 0);
                case 'CHRF1'
                    [~, ~, PARAM_OBJ{d}, HRF_OBJ{d}] = hrf_fit(data, HRF_PARAMS.TR, HRF_PARAMS.Condition{d}, HRF_PARAMS.T, HRF_PARAMS.types{t}, 1);
                case 'CHRF2'
                    [~, ~, PARAM_OBJ{d}, HRF_OBJ{d}] = hrf_fit(data, HRF_PARAMS.TR, HRF_PARAMS.Condition{d}, HRF_PARAMS.T, HRF_PARAMS.types{t}, 2);

                otherwise
                    error('No valid fit-type. Choose IL, FIR/sFIR or CHRF0/CHRF1/CHRF2')
            end

            
            for c=1:numel(HRF_PARAMS.Condition)
    
                HRF_OBJ{d}{c}.fullpath=sprintf([outfile, '_type-', HRF_PARAMS.types{t}, '_condition-', HRF_PARAMS.CondNames{c}, '_fit.nii']);
                PARAM_OBJ{d}{c}.fullpath=sprintf([outfile, '_type-', HRF_PARAMS.types{t}, '_condition-', HRF_PARAMS.CondNames{c}, '_params.nii']);
                try
                    write(HRF_OBJ{d}{c}, 'overwrite');
                    write(PARAM_OBJ{d}{c}, 'overwrite');
                catch
                    warning('Not able to write one or more files.');
                end
    
            end
    
            % HRF_local{d} = cell(1, numel(rois));
            % tc_local{d} = cell(1, numel(rois));

        end
    end

    % Generate an HRF and tc for every datafile and concatenate them
    % together.

    for d=1:numel(HRF_OBJ)
        for t=1:numel(HRF_OBJ{d})
            [HRF.fit{d,t}, tc{d,t}]=extractHRF(d{1}, [SPM.Sess(i).U.name], at, rois);
        end
    end

    % This results in an HRF struct: HRF(d, t, r, c), tc(d,t,r,c)
    

    delete(gcp('nocreate'));
end

function [tc, HRF, HRF_OBJ, PARAM_OBJ]=roiTS_fitHRF_SPM(SPM, HRF_PARAMS, rois, at, outfile)

    % Make directories for files if needed
    if ~isempty(fileparts(outfile))
        disp(fileparts(outfile));
        if ~exist(fileparts(outfile), 'dir')
            mkdir(fileparts(outfile));
        end
    end

    % Make sure CondNames are valid before continuing:
    HRF.atlas=at;
    HRF.region=rois;
    HRF.types=HRF_PARAMS.types;
    HRF.name=unique({SPM.xY.VY.fname}');

    if numel(HRF.name)>1
        for c=1:numel(HRF.name)
            % Initialize 'tc' and 'temp_HRF_fit' cell arrays
            tc{c} = cell(1, numel(HRF_PARAMS.types));
            temp_HRF_fit{c} = cell(1, numel(HRF_PARAMS.types));
        end
    else
        % Initialize 'tc' and 'temp_HRF_fit' cell arrays
        tc = cell(1, numel(HRF_PARAMS.types));
        temp_HRF_fit = cell(1, numel(HRF_PARAMS.types));
    end

    % Carve up SPM's design matrix for each image
        % DX=cell(1,numel(SPM.nscan));

    % HRF_PARAMS.CondNames=matlab.lang.makeValidName(HRF_PARAMS.CondNames);
    
    HRF_OBJ=cell(1,numel(HRF_PARAMS.types));
    PARAM_OBJ=cell(1,numel(HRF_PARAMS.types));
    info=cell(1,numel(HRF_PARAMS.types));
    
    % Initialize the parallel pool if it's not already running
    if isempty(gcp('nocreate'))
        parpool;
    end

    CondNames=cell(1, numel(HRF_PARAMS.types));

    parfor t=1:numel(HRF_PARAMS.types)
    % for t=1:numel(HRF_PARAMS.types)

    % First check to see if valid images have already been generated. Then
    % we don't have to regenerate them.
        for d = 1:numel(HRF.name)

            CondNames{t}{d} = strtrim([SPM.Sess(d).U.name]');
        
            for c=1:numel(CondNames{t}{d})
                HRF_OBJ_path=[outfile, '_type-', HRF_PARAMS.types{t}, '_condition-', CondNames{t}{d}{c}, '_fit.nii'];
                PARAM_OBJ_path=[outfile, '_type-', HRF_PARAMS.types{t}, '_condition-', CondNames{t}{d}{c}, '_params.nii'];
                if exist(HRF_OBJ_path, 'file')
                    HRF_OBJ{t}{d}{c}=fmri_data(HRF_OBJ_path);
                    PARAM_OBJ{t}{d}{c}=fmri_data(PARAM_OBJ_path);
                    files_exist=1;
                    % disp('found file')
                else
                    % disp('did not find file')
                    files_exist=0;
                end
            end
        end

        if ~files_exist
            % Generate the files otherwise.
            warning('off', 'all');
            if contains(HRF_PARAMS.types{t}, 'IL'), [~, ~, PARAM_OBJ{t}, HRF_OBJ{t}] = hrf_fit(SPM, HRF_PARAMS.T, 'IL', 0);, end
            if contains(HRF_PARAMS.types{t}, 'FIR'), [~, ~, PARAM_OBJ{t}, HRF_OBJ{t}, info{t}] = hrf_fit(SPM, HRF_PARAMS.T, 'FIR', 0);, end
            if contains(HRF_PARAMS.types{t}, 'sFIR'), [~, ~, PARAM_OBJ{t}, HRF_OBJ{t}, info{t}] = hrf_fit(SPM, HRF_PARAMS.T, 'FIR', 1);, end
            if contains(HRF_PARAMS.types{t}, 'CHRF0'), [~, ~, PARAM_OBJ{t}, HRF_OBJ{t}, info{t}] = hrf_fit(SPM, HRF_PARAMS.T, 'CHRF', 0);, end
            if contains(HRF_PARAMS.types{t},'CHRF1'), [~, ~, PARAM_OBJ{t}, HRF_OBJ{t}, info{t}] = hrf_fit(SPM, HRF_PARAMS.T, 'CHRF', 1);, end
            if contains(HRF_PARAMS.types{t}, 'CHRF2'), [~, ~, PARAM_OBJ{t}, HRF_OBJ{t}, info{t}] = hrf_fit(SPM, HRF_PARAMS.T, 'CHRF', 2);, end
    
            if ~ismember(HRF_PARAMS.types{t}, {'IL', 'FIR', 'sFIR','CHRF0','CHRF1','CHRF2'})
                error('Not a valid fit-type. Choose IL, FIR/sFIR or CHRF0/CHRF1/CHRF2')
            end
    
            % The issue with this one is that now HRF_OBJ{t} and PARAM_OBJ{t}
            % are cell arrays with cells for each data file d.
    
            for d = 1:numel(info{t})
    
                CondNames{t}{d} = strtrim(info{t}{d}.names);
            
                for c=1:numel(CondNames{t}{d})
        
                    HRF_OBJ{t}{d}{c}.fullpath=[outfile, '_type-', HRF_PARAMS.types{t}, '_condition-', CondNames{t}{d}{c}, '_fit.nii'];
                    PARAM_OBJ{t}{d}{c}.fullpath=[outfile, '_type-', HRF_PARAMS.types{t}, '_condition-', CondNames{t}{d}{c}, '_params.nii'];
                    try
                        write(HRF_OBJ{t}{d}{c}, 'overwrite');
                        write(PARAM_OBJ{t}{d}{c}, 'overwrite');
                    catch
                        warning('Not able to write one or more files.');
                    end
        
                end
    
            end
        end
    end

    HRF.CondNames=CondNames{1};

    % This will return an HRF_OBJ that is sectioned by:
    % HRF Type (t), Run(d), Region (r), Condition (c)

    % Rearrange HRF_OBJ{t, d, c} and PARAM_OBJ{t, d, c}



    % Generate an HRF and tc for every HRF_OBJ datafile (type x session x Condition) and concatenate them
    % together.

    HRF_arr=cell(1,numel(HRF_OBJ));
    for d=1:numel(HRF.name) % Organize by Session, then type, then condition:
        HRF_struct=HRF;
        HRF_struct.name=HRF.name{d};
        HRF_struct.CondNames=HRF.CondNames{d};
        tic
        for t=1:numel(HRF.types)
            [HRF_struct.fit{t}, tc{d}{t}]=extractHRF(HRF_OBJ{t}{d}, [SPM.Sess(d).U.name], at, rois);
        end
        toc

        HRF_arr{d}=HRF_struct;
    end

    HRF=HRF_arr;

     % This results in an HRF struct: HRF(d, t, r, c), tc(d,t,r,c)



end