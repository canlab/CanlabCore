function out = searchlight_disti_Lukas(dat, dist_i, additional_inputs)
%% ========================================================================
% Accelerated, Parfor-Safe Group-Level Searchlight Decoding
% - Index-based sphere caching (RAM-safe)
% - Subject-wise permutation inference (paired A/B flips)
% - Optional fixed C (from whole-brain nested CV)
% - fmri_data reconstruction fully corrected (no remove_empty errors)
% - Right-tailed voxelwise permutation p-values (acc & AUC)
%
% IMPORTANT:
%   dat MUST ALREADY be masked and trimmed before calling this function.
% ========================================================================

%% ------------------------------------------------------------------------
% Parse optional inputs
% -------------------------------------------------------------------------
do_cross = false;
r = 3;
P = 500;
fixedC = []; % NEW: optional user-supplied C

for k = 1:length(additional_inputs)
   if ischar(additional_inputs{k})
       switch additional_inputs{k}

           case 'dat2'
               do_cross = true;
               dat2 = additional_inputs{k+1};
               additional_inputs{k}=[]; additional_inputs{k+1}=[];

           case 'algorithm_name'
               algorithm_name = additional_inputs{k+1};
               additional_inputs{k}=[]; additional_inputs{k+1}=[];

           case 'r'
               r = additional_inputs{k+1};
               additional_inputs{k}=[]; additional_inputs{k+1}=[];

           case 'cv_assign'
               cv_assign = additional_inputs{k+1};
               additional_inputs{k}=[]; additional_inputs{k+1}=[];

           case 'n_permutations'
               P = additional_inputs{k+1};
               additional_inputs{k}=[]; additional_inputs{k+1}=[];

           case 'fixedC'  % NEW: fixed SVM C value
               fixedC = additional_inputs{k+1};
               additional_inputs{k}=[]; additional_inputs{k+1}=[];
       end
   end
end

%% ------------------------------------------------------------------------
% Convert to single precision (faster)
% -------------------------------------------------------------------------
dat.dat = single(dat.dat);
dat.volInfo.xyzlist = single(dat.volInfo.xyzlist);

if do_cross
   dat2.dat = single(dat2.dat);
   dat2.volInfo.xyzlist = single(dat2.volInfo.xyzlist);
end

%% ------------------------------------------------------------------------
% dist_i corresponds to ALREADY masked+trimmed dat
% -------------------------------------------------------------------------
vox_to_run = find(dist_i(:))';
N = numel(vox_to_run);
nvox = dat.volInfo.n_inmask;

%% ------------------------------------------------------------------------
% Parfor constants
% -------------------------------------------------------------------------
Cdat   = parallel.pool.Constant(dat);
Calg   = parallel.pool.Constant(algorithm_name);
Ccv    = parallel.pool.Constant(cv_assign);
CfixedC = parallel.pool.Constant(fixedC);  % NEW

if do_cross
   Cdat2 = parallel.pool.Constant(dat2);
else
   emptydat = dat;
   emptydat.dat = [];
   emptydat.removed_voxels = [];
   Cdat2 = parallel.pool.Constant(emptydat);
end

xyz = dat.volInfo.xyzlist;

%% ========================================================================
% SPHERE CACHING (index-based)
% ========================================================================
fprintf('Caching %d spheres (radius=%d)...\n', N, r);

sphere_cache = cell(N,1);
for ii = 1:N
   c = vox_to_run(ii);
   dx = xyz(:,1) - xyz(c,1);
   dy = xyz(:,2) - xyz(c,2);
   dz = xyz(:,3) - xyz(c,3);
   sphere_cache{ii} = int32(find((dx.*dx + dy.*dy + dz.*dz) <= r*r));
end
Cspheres = parallel.pool.Constant(sphere_cache);

%% ========================================================================
% TEMPLATE fmri_data OBJECT (prevents huge copying inside parfor)
% ========================================================================
template = dat;
template.dat = [];
template.removed_voxels = false(size(dat.removed_voxels));
template.volInfo.xyzlist = [];
template.volInfo.wh_inmask = [];
template.volInfo.n_inmask = [];
Ctemplate = parallel.pool.Constant(template);

%% ========================================================================
% STEP 1 — REAL SEARCHLIGHT DECODING
% ========================================================================
fprintf('Computing REAL searchlight (%d voxels)...\n', N);
q_real = makeProgressTracker(N);

real_acc = nan(N,1,'single');
real_auc = nan(N,1,'single');
real_se = nan(N,1,'single');
real_r  = nan(N,1,'single');

parfor idx = 1:N

   base = Cdat.Value;
   keep = Cspheres.Value{idx};

   % --- FAST sphere-specific fmri_data object ---
   dat_local = Ctemplate.Value;
   dat_local.dat = base.dat(keep,:);
   rv = true(size(base.removed_voxels));
   rv(keep) = false;
   dat_local.removed_voxels = rv;

   dat_local.volInfo.xyzlist = base.volInfo.xyzlist(keep,:);
   dat_local.volInfo.wh_inmask = find(~rv);
   dat_local.volInfo.n_inmask = numel(dat_local.volInfo.wh_inmask);

   R = predict_center_fast(dat_local, Calg.Value, Ccv.Value, CfixedC.Value);

   real_acc(idx) = R.acc;
   real_auc(idx) = R.AUC;
   real_se(idx) = R.se;
   real_r(idx)  = R.r;

   send(q_real,1);
end

%% ========================================================================
% STEP 2 — SUBJECT-WISE PERMUTATIONS
% ========================================================================
fprintf('Running %d permutations...\n', P);
q_perm = makeProgressTracker(P);

Nsub = length(Cdat.Value.Y)/2;

perm_count_acc = zeros(N,1,'uint32');
perm_count_auc = zeros(N,1,'uint32');

parfor p = 1:P

   base = Cdat.Value;
   dat_perm = base;

   % --- subject-wise paired sign-flips ---
   flipvec = rand(Nsub,1) > 0.5;
   for s = 1:Nsub
       a = s; b = s + Nsub;
       if flipvec(s)
           dat_perm.Y(a) = -1;
           dat_perm.Y(b) = 1;
           tmp = dat_perm.dat(a,:);
           dat_perm.dat(a,:) = dat_perm.dat(b,:);
           dat_perm.dat(b,:) = tmp;
       else
           dat_perm.Y(a) = 1;
           dat_perm.Y(b) = -1;
       end
   end

   local_acc = zeros(N,1,'uint32');
   local_auc = zeros(N,1,'uint32');

   for idx = 1:N

       keep = Cspheres.Value{idx};

       dat_local = Ctemplate.Value;
       dat_local.dat = dat_perm.dat(keep,:);
       rv = true(size(dat_perm.removed_voxels));
       rv(keep) = false;
       dat_local.removed_voxels = rv;

       dat_local.volInfo.xyzlist = base.volInfo.xyzlist(keep,:);
       dat_local.volInfo.wh_inmask = find(~rv);
       dat_local.volInfo.n_inmask = numel(dat_local.volInfo.wh_inmask);

       R = predict_center_fast(dat_local, Calg.Value, Ccv.Value, CfixedC.Value);

       if R.acc >= real_acc(idx), local_acc(idx) = 1; end
       if R.AUC >= real_auc(idx), local_auc(idx) = 1; end
   end

   perm_count_acc = perm_count_acc + local_acc;
   perm_count_auc = perm_count_auc + local_auc;

   send(q_perm,1);
end

p_acc = (perm_count_acc + 1)/(P+1);
p_auc = (perm_count_auc + 1)/(P+1);

%% ========================================================================
% STEP 3 — BUILD OUTPUT
% ========================================================================
out = struct();
out.test_results = cell(1,1);

full_acc = nan(nvox,1,'single');
full_auc = nan(nvox,1,'single');
full_se  = nan(nvox,1,'single');
full_r   = nan(nvox,1,'single');
full_pacc = nan(nvox,1,'single');
full_pauc = nan(nvox,1,'single');

full_acc(vox_to_run) = real_acc;
full_auc(vox_to_run) = real_auc;
full_se(vox_to_run)  = real_se;
full_r(vox_to_run)   = real_r;
full_pacc(vox_to_run) = p_acc;
full_pauc(vox_to_run) = p_auc;

out.test_results{1}.acc       = full_acc;
out.test_results{1}.AUC       = full_auc;
out.test_results{1}.se        = full_se;
out.test_results{1}.r         = full_r;
out.test_results{1}.p_acc_perm = full_pacc;
out.test_results{1}.p_auc_perm = full_pauc;

end

%% ========================================================================
% HELPER FUNCTIONS
% ========================================================================

function R = predict_center_fast(dat_local, alg, cv, fixedC)
   [test_Y, dat_local] = setup_testvar(dat_local);

   if ~isempty(fixedC)
       [~, stats] = predict(dat_local, ...
           'algorithm_name', alg, ...
           'C', fixedC, ...
           'nfolds', cv, ...
           'verbose', 0);
   else
       [~, stats] = predict(dat_local, ...
           'algorithm_name', alg, ...
           'nfolds', cv, ...
           'verbose', 0);
   end

   R = get_test_results(stats, test_Y, alg);
end

function R = get_test_results(stats, test_Y, alg)
   Y = test_Y{1};
   valid = (Y ~= 0);
   ytrue = (Y(valid)==1);

   if contains(alg,'svm')
       scores = stats.dist_from_hyperplane_xval(valid);
   else
       scores = stats.yfit(valid);
   end

   acc = mean((scores>0)==ytrue);

   try
       [~,~,~,AUC] = perfcurve(ytrue, scores, true);
   catch
       AUC = NaN;
   end

   try
       pos = scores(ytrue==1);
       neg = scores(ytrue==0);
       n1 = numel(pos); n0 = numel(neg);
       Q1 = AUC/(2-AUC);
       Q2 = 2*AUC*AUC/(1+AUC);
       varAUC = (AUC*(1-AUC)+(n1-1)*(Q1-AUC^2)+(n0-1)*(Q2-AUC^2))/(n1*n0);
       seAUC = sqrt(max(varAUC,0));
   catch
       seAUC = NaN;
   end

   R.acc = acc;
   R.AUC = AUC;
   R.se = seAUC;
   R.r  = NaN;
end

function [test_Y, dat] = setup_testvar(dat)
   test_Y{1} = dat.Y;
   dat.Y = dat.Y(:,1);
end

function q = makeProgressTracker(totalCount)
   progress.total = totalCount;
   progress.current = 0;
   progress.lastUpdate = tic;

   q = parallel.pool.DataQueue;
   afterEach(q,@(~)update());

   function update()
       progress.current = progress.current + 1;
       if toc(progress.lastUpdate) > 0.1 || progress.current == progress.total
           fprintf('\rProgress: %d/%d (%.1f%%)', ...
               progress.current, progress.total, ...
               100*progress.current/progress.total);
           progress.lastUpdate = tic;
           if progress.current == progress.total, fprintf('\n'); end
       end
   end
end