function [SPMinfo, dat] = canlab_glm_extract_trial_info_and_betas(input_path)
% Meta-data extraction and fmri_data creation for single-trial models
% Extracts information and single-trial data table from an SPM.mat file
%
% [SPMinfo, dat] = canlab_glm_extract_info_and_fmri_data(input_path)
%
% Inputs:
% -----------------------------------------------------------------------
% input_path: Path string specifying single-trial model with SPM.mat
%
% Outputs:
% -----------------------------------------------------------------------
% SPMinfo: Meta-data structure defining models
% dat: fmri_data object with extracted beta images for event-related regressors
%      - dat.additional_info has single trial data table
%
% Features:
% -----------------------------------------------------------------------
% - Extracts a space-efficient set of meta-data about the model
% - Extracts a table of info about the event-related regressors (trials)
% - Sorts betas into chronological order in the fmri_data object
% - Unzips SPM.mat and image files if necessary, and re-zips them afterwards
%
% - Currently especially set up for single trial models; could be extended to others, and contrasts
% - Does not yet handle contrasts, but easy to extend
% - Sorts into chrono order, but also retains the column numbers in the table (single_trial_info_table), so the original order is recoverable.
% - Currently re-zips SPM.mat and deletes original only if SPM.mat was
% zipped with gz.  Does not delete original image files after re-zipping
% yet - this feature could be added
%
%
% Tor Wager, 6/25/2017

% LOAD SPM FILE
% - make sure we have the correct one
% - unzip if needed

SPMinfo = [];
dat = [];

spmfilename = fullfile(input_path, 'SPM.mat');

waszipped = false;

if ~exist(spmfilename, 'file')
    
    spmzipped = fullfile(input_path, 'SPM.mat.gz');
    
    if exist(spmzipped, 'file')
        
        waszipped = true;
        disp('Unzipping SPM.mat.gz');
        gunzip(spmzipped)
        
    end
    
end

if exist(spmfilename, 'file')
    SPM = load(spmfilename);
    SPM = SPM.SPM;
    
    extracted_from_dir = fileparts(spmfilename);
    
    fprintf('Loaded %s\n', spmfilename);
    
    if waszipped
        disp('Re-zipping SPM.mat');
        gzip(spmfilename)
        delete(spmfilename)
    end
    
else
    disp('No SPM.mat file');
    return
end

%% Save minimal useful data in a lightweight structure

SPMinfo = struct('extracted_from_dir', extracted_from_dir);

% SPM fields to save
to_save = {'swd'                'xY.RT'      'xY.P'          'nscan'             'xBF' 'Sess' 'xX.X' 'xX.xKXs.X' 'xX.name'           'xM.VM.fname' 'xGX'                 'xVi.form'          'xCon'};

save_as = {'original_spm_dir'   'TR'         'input_images' 'images_per_session' 'xBF' 'Sess' 'X'    'Xfiltered' 'regressor_names'   'maskname' 'global_scaling_info' 'error_structure_form' 'xCon'};

for i = 1:length(to_save)
    
    spmfield = sprintf('SPM.%s', to_save{i});
    
    try
        
        mydata = eval(spmfield);
        SPMinfo.(save_as{i}) = mydata;
        
    catch
        fprintf('Field not found: %s\n', spmfield);
    end
    
end

%% Additional into

% Size of design matrix
[SPMinfo.n, SPMinfo.k] = size(SPMinfo.X);

SPMinfo.nsess = length(SPMinfo.images_per_session);

SPMinfo.hpfilter_by_session = cat(1, SPM.xX.K.HParam);

SPMinfo.TR_by_session = cat(1, SPM.xX.K.RT);

% List of betas
SPMinfo.beta_image_names = char(SPM.Vbeta.fname);


%% Get table of onsets, durations, names, beta images for single trials
% - For events in event-related design, with U field
% -----------------------------------------------------------------------

dat = []; % for fmri_data object


% Figure out if this is a single-trial model
% - for now, assume all events are single events.
% - needs code update to handle mixed single-trial and other events
% -----------------------------------------------------------------------
clear is_single_trial

k = SPMinfo.nsess;
is_single_trial = zeros(1, k);

for i = 1:k
    
    if ~isempty(SPM.Sess(i).U)
        % We have trial onsets
        is_single_trial(i) = length(cat(1, SPM.Sess(i).U.ons)) == length(SPM.Sess(i).U);
        
    else
        is_single_trial(i) = 1;  % default to 1 in case only some runs contain trials
        
    end
end

is_single_trial = all(is_single_trial);

%% Create single-trial table
% -----------------------------------------------------------------------

if is_single_trial
    
    clear ons_names ons durs sess_col continuous_regs continuous_reg_names wh_spike trial_num sess_trial_betas
    
    [sess_col, ons_names, ons, durs, trial_num, sess_trial_betas, sess_trial_image_names, run_num, continuous_regs, continuous_reg_names, wh_spike] = deal(cell(1, k));
    
    for i = 1:k
        % For each session
        % ---------------------------------------------------------------------
        
        sess_col{i} = SPM.Sess(i).col';  % columns in X for this session
        
        if ~isempty(SPM.Sess(i).U)
            % We have trial onsets
            
            ons_names{i} = cat(1, SPM.Sess(i).U.name);
            
            ons{i} = cat(1, SPM.Sess(i).U.ons);
            
            durs{i} = cat(1, SPM.Sess(i).U.dur);
            
            % trials are not in chronological order. Get which chronological
            % position each trial is:
            
            trial_num{i} = rankdata(ons{i});
            
            % Column number in xX for trials
            sess_trial_betas{i} = sess_col{i}(1:length(ons{i}));
            
            for j = 1:length(ons{i})
                
                sess_trial_image_names{i}{j, 1} = sprintf('beta_%04d.img', sess_trial_betas{i}(j));
                
            end
            
            % run number (for data table)
            
            run_num{i} = repmat(i, length(ons{i}), 1);
            
        end
        
        % Continuous regressors
        continuous_regs{i} = SPM.Sess(i).C.C;
        
        continuous_reg_names{i} = SPM.Sess(i).C.name;
        
        % Single trials: Which have exactly one onset of value 1
        wh_spike{i} =-all(continuous_regs{i} == 1 | continuous_regs{i} == 0, 1) & sum(continuous_regs{i}) == 1;
        
    end % Session
    
    % Concatenate trial info
    
    single_trial_info_table = table();
    single_trial_info_table.image_names = cat(1, sess_trial_image_names{:});
    
    single_trial_info_table.condition = cat(1, ons_names{:});
    
    single_trial_info_table.ons = cat(1, ons{:});
    single_trial_info_table.durs = cat(1, durs{:});
    single_trial_info_table.trial_num = cat(1, trial_num{:});
    single_trial_info_table.beta_num = cat(1, sess_trial_betas{:});
    single_trial_info_table.run_num = cat(1, run_num{:});
    
    % Sort into chronological order
    
    single_trial_info_table = sortrows(single_trial_info_table, {'run_num' 'trial_num'}, 'ascend');
    
    
    SPMinfo.single_trial_info_table = single_trial_info_table;
    
    
    %% Extract data from images
    
    cd(extracted_from_dir)
    
    % Unzip
    !gunzip -f *.img.gz
    !gunzip -f *.hdr.gz
    !gunzip -f *.nii.gz
    
    dat = fmri_data(single_trial_info_table.image_names);
    dat = enforce_variable_types(dat);
    
    % Store meta-data with object in additional_info field:
    dat.additional_info = single_trial_info_table;
    
    % zip - do later outside this fcn
%     !gzip *.img
%     !gzip *.nii
    
    
end % if is single trial model



end % function






%%
%
%
%
% %% Get single trial betas
%
% clear ons_names ons durs sess_col continuous_regs continuous_reg_names wh_spike trial_num sess_trial_betas
% k = length(SPM.Sess);
%
% for i = 1:k
%     % For each session
%     % ---------------------------------------------------------------------
%
%     sess_col{i} = SPM.Sess(i).col';
%
%     if ~isempty(SPM.Sess(i).U)
%         % We have trial onsets
%
%         ons_names{i} = cat(1, SPM.Sess(i).U.name);
%
%         ons{i} = cat(1, SPM.Sess(i).U.ons);
%
%         durs{i} = cat(1, SPM.Sess(i).U.dur);
%
%         % trials are not in chronological order. Get which chronological
%         % position each trial is:
%
%         trial_num{i} = rankdata(ons{i});
%
%         % Column number in xX for trials
%         sess_trial_betas{i} = sess_col{i}(1:length(ons{i}));
%
%         for j = 1:length(ons{i})
%
%             sess_trial_image_names{i}{j, 1} = sprintf('beta_%04d.img', sess_trial_betas{i}(j));
%
%         end
%
%         % run number (for data table)
%
%         run_num{i} = repmat(i, length(ons{i}), 1);
%
%     end
%
%     % Continuous regressors
%     continuous_regs{i} = SPM.Sess(i).C.C;
%
%     continuous_reg_names{i} = SPM.Sess(i).C.name;
%
%     % Single trials: Which have exactly one onset of value 1
%     wh_spike{i} =-all(continuous_regs{i} == 1 | continuous_regs{i} == 0, 1) & sum(continuous_regs{i}) == 1;
%
% end % Session
%
% % Concatenate trial info
%
% single_trial_info_table = table();
% single_trial_info_table.image_names = cat(1, sess_trial_image_names{:});
%
% single_trial_info_table.condition = cat(1, ons_names{:});
%
% single_trial_info_table.ons = cat(1, ons{:});
% single_trial_info_table.durs = cat(1, durs{:});
% single_trial_info_table.trial_num = cat(1, trial_num{:});
% single_trial_info_table.beta_num = cat(1, sess_trial_betas{:});
% single_trial_info_table.run_num = cat(1, run_num{:});
%
% % Sort into chronological order
%
% single_trial_info_table = sortrows(single_trial_info_table, {'run_num' 'trial_num'}, 'ascend');
%
% %% Extract data from images
%
% dat = fmri_data(single_trial_info_table.image_names);
% dat = enforce_variable_types(dat);
%
% % Store meta-data with object in additional_info field:
% dat.additional_info = single_trial_info_table;
%
%
%
% %% Sizes of fields in SPM.mat
%
% % Most of the space in SPM.mat files in in SPM.xY.VY(1).private
% % Other big fields: SPM.xVi.CY;  SPM.xX.V;  SPM.xX.xKXs; SPM.xX.xKXs.u;
% % SPM.xX.xKXs.X;
%
% % n = fieldnames(SPM);
% % for i = 1:length(n)
% % x = SPM.(n{i});
% % ww{i} = whos('x');
% % end
%
% % xY	216667520
% % xBF	4798
% % nscan	 56
% % Sess	1669204
% % xGX	50672
% % xX	79946124
% % xVi	76756762
% % SPMid	 42
% % xM	62815
% % xsDes	1576
% % swd	284
% % xVol	5151999
% % Vbeta	5042010
% % VResMS	11999
% % VM	11995
% % xCon	  0