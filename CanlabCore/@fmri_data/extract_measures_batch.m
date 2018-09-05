function DAT = extract_measures_batch(data_obj)

t1 = clock;
dashes = '----------------------------------------------';
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);


%% Initialize DAT output
% data_obj = DATA_OBJ{1};

% Save some meta-data about the images

DAT = [];
DAT.extracted_on_date = scn_get_datetime;
DAT.image_names = data_obj.image_names;
DAT.fullpath = data_obj.fullpath;

%% Mahal

printhdr('Mahalanobis distance for each image, based on correlation');

[mahal_dist, expected_dist, p_val, wh_outlier_uncorr, wh_outlier_corr] = mahal(data_obj, 'noplot');

DAT.mahalanobis = table(mahal_dist, expected_dist, p_val, wh_outlier_uncorr, wh_outlier_corr);

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
    
    DAT.PARCELS.(parcellation_name).conditions = DAT.conditions;
    
    k = length(DAT.conditions);
    
    for i = 1:k
        
        % Extract mean values
        
        parcel_means = apply_parcellation(DATA_OBJ{i}, parcel_obj);
        
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
        
        printstr(signame)
        
        for i = 1:k
            
            % Extract mean values
            
            local_pattern = apply_parcellation(DATA_OBJ{i}, parcel_obj, 'pattern_expression', sig_obj);
            
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

