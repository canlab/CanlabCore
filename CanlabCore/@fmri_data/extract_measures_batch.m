function DAT = extract_measures_batch(data_obj)
% Extracts a set of measures relevant for pattern-based and network-based analyses
%
% DAT = extract_measures_batch(data_obj)
% 
% This is a method called extract_measures_batch for fmri_objects.  
% It Extracts a set of measures relevant for pattern-based and network-based analyses. 
% The idea is to aggregate these across studies, and pull relevant measures from the set 
% for particular analyses.  It returns the following structure:
%
% Tor Wager, August 2018
%
% 
% DAT = 
% 
%   struct with fields:
% 
%        extracted_on_date: '05-Sep-2018_02_49? 	Date information was extracted
%              image_names: 'wrrest_mb8_r1.nii? 	Original image names (no paths)
%                 fullpath: [914×119 char]          Full path names for all volumes (for provenance)
%              mahalanobis: [914×5 table]           Mahalanobis distances (for weighting/nuisance)and outlier ID logical
%                    rmssd: [1×1 struct] 			Root mean square successive differences and outlier ID logical
%     gray_white_csf_table: [914×5 table]           Global avg gray, white, CSF, and 5 principal components	 for each; for nuisance
%                  npsplus: [1×1 struct] 			Multivariate pattern responses for CANlab measures (NPS, more)
%            kragelemotion: [1×1 struct]            Multivariate pattern responses for Kragel 2015 emotion classification
%                 kragel18: [1×1 struct] 			Multivariate PLS pattern responses for Kragel 2018 Nat Neurosci and subregions
%                 pain_pdm: [1×1 struct] 			Multivariate pattern responses for Geuter et al.?s Prin. Dirs of Medation (10 patterns, and combined)
%                  PARCELS: [1×1 struct] 			Parcellations: Averages for each parcel, and local pattern responses (selected)
%
% Tables:
% Some variables are in Matlab table objects, e.g., DAT.mahalanobis. 
% 
% Access names like this:
% DAT.mahalanobis.Properties.VariableNames
% See methods(table) for more operations.  table2array() is useful for extracting data matrices.
% outliers = DAT.mahalanobis.wh_outlier_uncorr;  % returns a logical vector of outlier images id'd as having high Mahalanobis distances
%
% OUTLIERS
% To get a reasonable set of outlier image indicators (indicator vector):
% outliers = DAT.mahalanobis.wh_outlier_uncorr | DAT.rmssd.wh_outliers_rmssd;  % returns a logical vector of outlier images id'd as having high Mahalanobis distances
%
% EXTRACTED GLOBAL TISSUE COMPARTMENTS
% DAT.gray_white_csf_table.Properties.VariableNames
% 
% {'gwcsf'       }      % n_images x 3, global mean gray, white, CSF (assumes MNI space; group template)
% {'gray5'       }      % n_images x 5, first 5 principal components across gray matter
% {'white5'      }      % n_images x 5, first 5 principal components across white matter
% {'csf5'        }      % n_images x 5, first 5 principal components across CSF
% {'gwcsf_l2norm'}      % n_images x 3, L2 norms for each image x [gray white CSF]
%
% 
% NUISANCE COVARIATES:
% -------------------------------------------------------------------------
% Here is a reasonable set of nuisance covariates. If you are contenating
% across runs (when scaner starts and stops), add indicator vectors for
% run, plus movement covariates (e.g., 24 per run, not included.)
%
% First, turn outlier vector into a set of separate  dummy regressors for each outlier point
% outliers = DAT.mahalanobis.wh_outlier_uncorr | DAT.rmssd.wh_outliers_rmssd;  % returns a logical vector of outlier images id'd as having high Mahalanobis distances
% outliers = double(outliers);
% outliers(outliers > 0) = find(outliers);
% [outlier_indic, outlier_levels] = condf2indic(outliers);
% outlier_indic(:, outlier_levels == 0) = [];
% outlier_names = {};
% for i = 1:size(outlier_indic, 2)
%     outlier_names{i} = sprintf('Spike%3.0f', i);
% end
%
% Second, pull out other covariates and concatenate them into a matrix:
% cov_matrix = [DAT.gray_white_csf_table.csf5 DAT.gray_white_csf_table.gwcsf_l2norm(:, 3) zscore(DAT.rmssd.rmssd) outlier_indic];
% cov_names = [{'CSF_comp1' 'CSF_comp1' 'CSF_comp1' 'CSF_comp1' 'CSF_comp1' 'CSF_l2norm' 'rmssd'} outlier_names];
%
% SIGNATURES
% -------------------------------------------------------------------------
% A number of "signature patterns" are extracted, with fields indicating
% whether the data were scaled first or "raw" (preprocessed), and for
% different similarity metrics. 
% 
% DAT.npsplus.raw .     % raw means we have not scaled the images, e.g., divided by L2 norm
% 
% ans = 
% 
%   struct with fields:
% 
%      dotproduct: [1×1 struct] % These are different similarity metrics
%      cosine_sim: [1×1 struct]
%     correlation: [1×1 struct]
%
%     similarity_metric: 'dotproduct'
%         image_scaling: 'none'
%        signaturenames: {1×14 cell}
%        conditionnames: {'C__1'}
%                   NPS: [914×1 table]  % Each of these is a table object with a different signature
%                NPSpos: [914×1 table]
%                NPSneg: [914×1 table]
%                 SIIPS: [914×1 table]
%                 PINES: [914×1 table]
%             Rejection: [914×1 table]
%                   VPS: [914×1 table]
%           VPS_nooccip: [914×1 table]
%                   GSR: [914×1 table]
%                 Heart: [914×1 table]
%          FM_Multisens: [914×1 table]
%               FM_pain: [914×1 table]
%         Empathic_Dist: [914×1 table]
%         Empathic_Care: [914×1 table]
%         
% Access them and build a matrix like this:
% 
% mynames = DAT.npsplus.raw.dotproduct.signaturenames;
% my_signature_data = [];
% for i = 1:length(mynames)
%     my_signature_data(:, i) = table2array(DAT.npsplus.raw.dotproduct.(mynames{i}));
% end
%
% PARCELS
% -------------------------------------------------------------------------
%
% DAT.PARCELS . contains extracted data for two different parcellations.
% For both, it extracts both region averages from each parcel and local
% pattern expression from selected parcels * patterns.
% 
% ans = 
% 
%   struct with fields:
% 
%     canlab2018_2mm: [1×1 struct]  % ~500-region atlas composite from multiple published atlases and named ROIs
%      yeo17networks: [1×1 struct]  % 16 unique rsfMRI networks, separated into left and right hemispheres
%      
% DAT.PARCELS.yeo17networks
% 
% ans = 
% 
%   struct with fields:
% 
%        parcel_obj: [1×1 atlas]    % Original atlas object, with region labels, etc.
%             means: [1×1 struct]   % .dat has a time x parcels matrix of mean data from each parcel
%               NPS: [1×1 struct]   % Local pattern expression for the NPS in each parcel
%             SIIPS: [1×1 struct]
%             PINES: [1×1 struct]
%               VPS: [1×1 struct]
%     Empathic_Care: [1×1 struct]     
%
% EXAMPLE
% -------------------------------------------------------------------------
% % Load a sample set of test images
% test_images = load_image_set('emotionreg');           % Load some test images
%
% % Extract parcel data, signature responses, etc.
% DAT = extract_measures_batch(test_images);
%
% % Visualize parcels and the data in the whole-brain parcel matrix
% orthviews(DAT.PARCELS.canlab2018_2mm.parcel_obj);
% disp(DAT.PARCELS.canlab2018_2mm.parcel_obj.labels);
%
% create_figure('parcel_data'); 
% imagesc(DAT.PARCELS.canlab2018_2mm.means.dat{1});
% axis tight; ylabel('Image'); xlabel('Parcel number'); colorbar;


t1 = clock;
dashes = '----------------------------------------------';
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);


%% Initialize DAT output
% data_obj = DATA_OBJ{1};

% Save some meta-data about the images

DAT = [];
DAT.extracted_with = 'Extracted parcels, signatures, and local signatures with fmri_data.extract_measures_batch()';
DAT.extracted_on_date = scn_get_datetime;
DAT.image_names = data_obj.image_names;
DAT.fullpath = data_obj.fullpath;

%% Mahal

printhdr('Mahalanobis distance for each image, based on correlation');

[mahal_dist, expected_dist, p_val, wh_outlier_uncorr, wh_outlier_corr] = mahal(data_obj, 'noplot');

DAT.mahalanobis = table(mahal_dist, expected_dist, p_val, wh_outlier_uncorr, wh_outlier_corr);

%% RMSSD / dVARS

% root mean square successive differences (dvars)
% mydiffs = mean(diff(data_obj.dat') .^ 2, 2) .^ .5;
% DAT.rmssd = [0; mydiffs];

DAT.rmssd.descrip = 'Root mean square successive diffs in images. Outliers IDd at 3 sd above mean';

[~, DAT.rmssd.rmssd, DAT.rmssd.wh_outliers_rmssd] = preprocess(data_obj, 'outliers_rmssd');


%% Extract global gray, white, CSF

printhdr('Extracting global gray, white, CSF');

[gwcsf, components, ~, gwcsf_l2norm] = extract_gray_white_csf(data_obj);
[gray5, white5, csf5] = deal(components{:});

DAT.gray_white_csf_table = table(gwcsf, gray5, white5, csf5, gwcsf_l2norm);

%% Extract whole-brain signature patterns

printhdr('Extracting whole-brain signature patterns');

% SIG = apply_all_signatures(data_obj, ['similarity_metric', sim_metric, 'image_scaling', imgscaling, 'conditionnames', conditionnames, 'image_set', image_set_name);

% Dot product metric
DAT.npsplus.raw.dotproduct = apply_all_signatures({data_obj});
DAT.kragelemotion.raw.dotproduct = apply_all_signatures({data_obj}, 'image_set', 'kragelemotion');
DAT.kragel18.raw.dotproduct = apply_all_signatures({data_obj}, 'image_set', 'kragel18');
DAT.pain_pdm.raw.dotproduct = apply_all_signatures({data_obj}, 'image_set', 'pain_pdm');

% Cosine similarity
DAT.npsplus.raw.cosine_sim = apply_all_signatures({data_obj}, 'similarity_metric', 'cosine_similarity');
DAT.kragelemotion.raw.cosine_sim = apply_all_signatures({data_obj}, 'image_set', 'kragelemotion', 'similarity_metric', 'cosine_similarity');
DAT.kragel18.raw.cosine_sim = apply_all_signatures({data_obj}, 'image_set', 'kragel18', 'similarity_metric', 'cosine_similarity');
DAT.pain_pdm.raw.cosine_sim = apply_all_signatures({data_obj}, 'image_set', 'pain_pdm', 'similarity_metric', 'cosine_similarity');

% Correlation
DAT.npsplus.raw.correlation = apply_all_signatures({data_obj}, 'similarity_metric', 'correlation');
DAT.kragelemotion.raw.correlation = apply_all_signatures({data_obj}, 'image_set', 'kragelemotion', 'similarity_metric', 'correlation');
DAT.kragel18.raw.correlation = apply_all_signatures({data_obj}, 'image_set', 'kragel18', 'similarity_metric', 'correlation');
DAT.pain_pdm.raw.correlation = apply_all_signatures({data_obj}, 'image_set', 'pain_pdm', 'similarity_metric', 'correlation');

%% THIS section extracts region averages and local pattern expression in each
% parcel of the canlab2018_2mm atlas
%
% - extracts local parcel means and all local pattern expression for all signatures.
% - Preps file Parcellation_data.mat saved in "results" folder.
% - Allow us to compare local pattern expression across parcels and compare
%   local patterns within parcels.
%
% Saves variable named whatever the var "parcellation_name" is set to
% - Need to break files up into one file per parcellation due to memory space

% Get parcels
% --------------------------------------------------------------------

parcellation_names = {'canlab2018_2mm' 'yeo17networks'};

% Prep: load sig objects for local pattern expression
% --------------------------------------------------------------------

[signature_obj, signames] = load_image_set('npsplus', 'noverbose');

% Select a few signatures (save time extracting)
wh_sigs = [1 4 5 7 14];
signames = signames(wh_sigs);
signature_obj = get_wh_image(signature_obj, wh_sigs);


% --------------------------------------------
% Resample data to standard space 
% for compat. with atlases or sigs 
% Can save time later with resampling...
% --------------------------------------------

disp('Resampling data to space of signatures')
data_obj_sig_space = resample_space(data_obj, signature_obj);

ref_atlas_obj = load_atlas(parcellation_names{1});

disp('Resampling data to space of atlas object')
data_obj_atlas_space = resample_space(data_obj, ref_atlas_obj);

clear ref_atlas_obj

% isdiff = compare_space(data_obj, signature_obj);
% 
% if isdiff == 0 || isdiff == 3 % Same space, diff voxels
%    % do nothing
% elseif isdiff == 2
%     error('Invalid object: volInfo structure missing for one or more objects');
% else
%     data_obj_sig_space = resample_space(data_obj, signature_obj);
% end


for pn = 1:length(parcellation_names)
    
    
    tic
    disp('Loading parcels');
    
    parcellation_name = parcellation_names{pn};   %'canlab2018_2mm';
    parcel_obj = load_atlas(parcellation_name);
    
    printhdr(parcellation_name)
    
    % Space-saving
    parcel_obj = enforce_variable_types(parcel_obj);
    
    DAT.PARCELS.(parcellation_name) = [];
    DAT.PARCELS.(parcellation_name).parcel_obj = parcel_obj;
    % DAT.PARCELS.(parcellation_name).regions = region(parcel_obj, 'unique_mask_values');
    
    toc
    
    %% Get mean activation for each parcel
    % --------------------------------------------------------------------
    tic
    printhdr('Extracting parcel mean values');
    
    %DAT.PARCELS.(parcellation_name).conditions = DAT.conditions;
    
    k = 1; %length(DAT.conditions);
    
    for i = 1:k
        
        % Extract mean values
        
        parcel_means = apply_parcellation(data_obj_atlas_space, parcel_obj);
        
        DAT.PARCELS.(parcellation_name).means.dat{i} = parcel_means;
        
        % Do a t-test on each parcel
        [h, p, ci, stat] = ttest(double(parcel_means));
        
        DAT.PARCELS.(parcellation_name).means.group_t{i} = stat.tstat;
        DAT.PARCELS.(parcellation_name).means.group_p{i} = p;
        
    end
    
    % FDR-correct across all conditions and parcels
    
    all_p = cat(2, DAT.PARCELS.(parcellation_name).means.group_p{:});
    DAT.PARCELS.(parcellation_name).means.fdr_p_thresh = FDR(all_p, .05);
    
    for i = 1:k
        
        DAT.PARCELS.(parcellation_name).means.fdr_sig{i} = DAT.PARCELS.(parcellation_name).means.group_p{i} < DAT.PARCELS.(parcellation_name).means.fdr_p_thresh;
        
    end
    
    toc
    
    %% Get pattern response for each parcel
    % --------------------------------------------------------------------
    
    tic
    printhdr('Extracting parcel local signature patterns');
    
    for mysig = 1:length(signames)
        
        sig_obj = get_wh_image(signature_obj, mysig);
        
        signame = signames{mysig};
        signame = strrep(signame, '-', '_');
        signame = strrep(signame, ' ', '_');
        
        disp(signame)
        
        for i = 1:k
            
            % Extract mean values
            
            local_pattern = apply_parcellation(data_obj_atlas_space, parcel_obj, 'pattern_expression', sig_obj);
            
            DAT.PARCELS.(parcellation_name).(signame).dat{i} = local_pattern;
            
            % Do a t-test on each parcel
            [h, p, ci, stat] = ttest(double(local_pattern));
            
            DAT.PARCELS.(parcellation_name).(signame).group_t{i} = stat.tstat;
            DAT.PARCELS.(parcellation_name).(signame).group_p{i} = p;
            
        end
        
        % FDR-correct across all conditions and parcels
        
        all_p = cat(2, DAT.PARCELS.(parcellation_name).(signame).group_p{:});
        
        pthr = FDR(all_p, .05);
        if isempty(pthr), pthr = eps; end
        
        DAT.PARCELS.(parcellation_name).(signame).fdr_p_thresh = pthr;
        
        for i = 1:k
            
            DAT.PARCELS.(parcellation_name).(signame).fdr_sig{i} = DAT.PARCELS.(parcellation_name).(signame).group_p{i} < DAT.PARCELS.(parcellation_name).(signame).fdr_p_thresh;
            
        end
        
        toc
        
    end  % signature
    
    
end % parcellation

disp('Total elapsed time:')
etime(clock, t1)

end % function


%% Format into statistic_image objects for display, etc.
% ADD these to PARCEL structure
%
% display with, e.g., orthviews(DAT.PARCELS.(parcellation_name).(signame).t_statistic_obj{3})

%parcel_obj = DAT.PARCELS.(parcellation_name).parcel_obj;
%
% printhdr('Reconstructing parcel-wise t-statistic objects');
%
% for mysig = 1:length(signames)
%
%     signame = signames{mysig};
%     signame = strrep(signame, '-', '_');
%     signame = strrep(signame, ' ', '_');
%
%     printstr(signame)
%
%     DAT.PARCELS.(parcellation_name).(signame) = plugin_get_parcelwise_statistic_images(parcel_obj, DAT.PARCELS.(parcellation_name).(signame) );
%
% end
%
%

