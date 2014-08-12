
% needs output from classify_bayes
% indx = voxels (all voxels) x classes (task states)
% ptask = p(task) within reduced feature set

nclasses = size(ptask, 1);
whsave = sum(indx > 0, 2) > 0;

% Get MDS coordinates
% -----------------------------------------------
r = corrcoef(ptask');
D = 1 - r;

[mds_dims, stress, disp] = shepardplot(D, 4);


% Get coordinates of ideal post. probabilities
% -----------------------------------------------
prjn_matrix = pinv(mds_dims);

ideal_coords = prjn_matrix * eye(nclasses);


%% Plot them
% -----------------------------------------------
colors = {[1 0 0] [0 1 0] [1 0 1] [1 1 0] [0 0 1]};
%plot(mds_dims(:,1), mds_dims(:,2), 'ko', 'MarkerFaceColor', [.5 .5 .5]);
nmdsfig(ideal_coords', 'classes', 1:nclasses, 'names', names, 'colors', colors, 'sizes', 16);


%% Set up to classify images

meta_space_image = MC_Setup.V.fname;
% /Users/tor/Documents/matlab_code/3DheadUtility/canonical_brains/SPM99_sca
% lped_brains/scalped_avg152T1_graymatter_smoothed.img

all_images = spm_get(Inf);
nimgs = size(all_images, 1);


%% Take a new data vector and plot it in the space
% -----------------------------------------------
hold on;

for i = 1:nimgs
    image_to_classify = all_images(i, :);

    % /Users/tor/Documents/Tor_Documents/CurrentExperiments/Lab_Emotion/Resil2_Speech_Task/SPEECH_TASK_RESULTS/speech_seed_sh5_norob/avg_r_HR_0001.img
    imgdata = scn_map_image(image_to_classify, meta_space_image);
    imgvec = imgdata(MC_Setup.volInfo.wh_inmask);

    imgvec = imgvec(whsave);

    imgvec = scale(imgvec, 1);
    imgvec(isnan(imgvec)) = 0;

    [predicted_class, ml, lr, taskprob] = classify_choose_most_likely(ptask, imgvec);
    coords = prjn_matrix * taskprob;
    plot3(coords(1), coords(2), coords(3), 'o', 'Color', colors{predicted_class}, 'MarkerFaceColor', colors{predicted_class});

    fprintf(1, '%3.0f  ', predicted_class);
    fprintf(1, '%3.2f  ', taskprob);
    fprintf(1, '\n')
    
    drawnow
    pause(.5)

end