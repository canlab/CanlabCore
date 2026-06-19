function out = searchlight_disti_Lukas(dat, dist_i, additional_inputs)
%% ========================================================================
% searchlight_disti_Lukas_tfce
%
% Group-level searchlight decoding with TFCE-only permutation inference
% using a LOCAL implementation of classic TFCE (Smith & Nichols, 2009).
%
% No dependency on SPM / TDT / CoSMoMVPA TFCE implementations.
%
% ASSUMPTIONS:
%  - dat is already masked & trimmed
%  - exactly 2 images per subject: [A1..AN B1..BN]
%  - dat.Y = [1..1 -1..-1]
%
% OUTPUT:
%  out.real_auc
%  out.real_t
%  out.TFCE_real
%  out.TFCE_real_max
%  out.TFCE_null_max
%  out.p_TFCE_global
%  out.p_TFCE_voxel
%  out.tfce_stat_image
% ========================================================================

%% -------------------- Defaults --------------------
r          = 3;
P          = 100;
alg        = 'cv_svm';
fixedC     = [];
cv_assign  = [];

% TFCE defaults (Smith & Nichols)
tfce_H     = 2;
tfce_E     = 0.5;
tfce_conn  = 26;

%% -------------------- Parse inputs --------------------
for k = 1:length(additional_inputs)
   if ischar(additional_inputs{k})
       switch lower(additional_inputs{k})
           case 'r'
               r = additional_inputs{k+1};
           case 'algorithm_name'
               alg = additional_inputs{k+1};
           case 'fixedc'
               fixedC = additional_inputs{k+1};
           case 'cv_assign'
               cv_assign = additional_inputs{k+1};
           case 'n_permutations'
               P = additional_inputs{k+1};
           case 'tfce_h'
               tfce_H = additional_inputs{k+1};
           case 'tfce_e'
               tfce_E = additional_inputs{k+1};
           case 'tfce_conn'
               tfce_conn = additional_inputs{k+1};
       end
   end
end

if ~isempty(cv_assign)
   Ccv = parallel.pool.Constant(cv_assign);
else
   Ccv = [];
end

%% -------------------- Dimensions --------------------
Y     = dat.Y;
Nsub  = numel(Y) / 2;
xyz   = single(dat.volInfo.xyzlist);

vox_to_run = find(dist_i(:));
Nvox = numel(vox_to_run);

fprintf('Searchlight voxels: %d\n', Nvox);
fprintf('Subjects: %d\n', Nsub);
fprintf('TFCE permutations: %d\n', P);

%% -------------------- Sphere caching --------------------
fprintf('Caching spheres (r=%d)...\n', r);
sphere_cache = cell(Nvox,1);
for v = 1:Nvox
   c = vox_to_run(v);
   dx = xyz(:,1) - xyz(c,1);
   dy = xyz(:,2) - xyz(c,2);
   dz = xyz(:,3) - xyz(c,3);
   sphere_cache{v} = int32(find(dx.^2 + dy.^2 + dz.^2 <= r^2));
end
Cspheres = parallel.pool.Constant(sphere_cache);

%% -------------------- fmri_data template --------------------
template = dat;
template.dat = [];
template.removed_voxels = false(size(dat.removed_voxels));
Ctemplate = parallel.pool.Constant(template);

%% ========================================================================
% STEP 1 — REAL SEARCHLIGHT
%% ========================================================================
fprintf('Running REAL searchlight decoding...\n');

tracker_real = ProgressTracker(Nvox);
q_real = parallel.pool.DataQueue;
afterEach(q_real,@(~)tracker_real.update());

real_auc = nan(Nvox,1,'single');

parfor v = 1:Nvox
   keep = Cspheres.Value{v};
   base = dat;

   dloc = Ctemplate.Value;
   dloc.dat = base.dat(keep,:);
   rv = true(size(base.removed_voxels));
   rv(keep)=false;
   dloc.removed_voxels = rv;

   if isempty(fixedC)
       if isempty(Ccv)
           [~,stats] = predict(dloc,'algorithm_name',alg,'verbose',0);
       else
           [~,stats] = predict(dloc,'algorithm_name',alg,'cv_assign',Ccv.Value,'verbose',0);
       end
   else
       if isempty(Ccv)
           [~,stats] = predict(dloc,'algorithm_name',alg,'C',fixedC,'verbose',0);
       else
           [~,stats] = predict(dloc,'algorithm_name',alg,'C',fixedC,'cv_assign',Ccv.Value,'verbose',0);
       end
   end

   scores = stats.dist_from_hyperplane_xval;
   ytrue = (Y==1);
   [~,~,~,AUC] = perfcurve(ytrue,scores,true);

   real_auc(v)=AUC;
   send(q_real,1);
end

%% -------------------- AUC -> t --------------------
real_t = auc_to_t(real_auc,Nsub,Nsub);

%% -------------------- TFCE on real map --------------------
vol_dim = dat.volInfo.dim;
tvol = zeros(vol_dim,'single');
tvol(~dat.removed_voxels)=real_t;

TFCE_real_vol = tfce_transform_3d(tvol,tfce_H,tfce_E,tfce_conn);
TFCE_real = TFCE_real_vol(~dat.removed_voxels);
TFCE_real_max = max(TFCE_real);

%% ========================================================================
% STEP 2 — TFCE-ONLY PERMUTATIONS
%% ========================================================================
fprintf('Running TFCE-only permutations...\n');

tracker_perm = ProgressTracker(P);
q_perm = parallel.pool.DataQueue;
afterEach(q_perm,@(~)tracker_perm.update());

maxTFCE_null = nan(P,1,'single');

parfor p=1:P
   perm_dat = dat;
   flip = rand(Nsub,1)>0.5;
   for s=1:Nsub
       a=s; b=s+Nsub;
       if flip(s), perm_dat.Y([a b])=perm_dat.Y([b a]); end
   end

   perm_auc = nan(Nvox,1,'single');
   for v=1:Nvox
       keep = Cspheres.Value{v};
       dloc = Ctemplate.Value;
       dloc.dat = perm_dat.dat(keep,:);
       perm_auc(v)=get_auc(Ctemplate.Value,perm_dat,keep,alg,fixedC,Ccv);
   end

   t_perm = auc_to_t(perm_auc,Nsub,Nsub);
   tvol_perm = zeros(vol_dim,'single');
   tvol_perm(~dat.removed_voxels)=t_perm;

   TFCE_perm = tfce_transform_3d(tvol_perm,tfce_H,tfce_E,tfce_conn);
   maxTFCE_null(p)=max(TFCE_perm(~dat.removed_voxels));
   send(q_perm,1);
end

%% -------------------- TFCE p-maps --------------------
p_TFCE_global = (sum(maxTFCE_null>=TFCE_real_max)+1)/(P+1);

TFCE_p_voxel = nan(size(TFCE_real),'single');
for v=1:numel(TFCE_real)
   TFCE_p_voxel(v)=(sum(maxTFCE_null>=TFCE_real(v))+1)/(P+1);
end

tfce_stat_img = statistic_image('type','p');
tfce_stat_img.volInfo = dat.volInfo;
tfce_stat_img.removed_voxels = dat.removed_voxels;
tfce_stat_img.dat = TFCE_p_voxel;
tfce_stat_img.dat_descrip = ['uncorrected TFCE p-values based on ' num2str(P) ' permutations'];
tfce_stat_img.p = TFCE_p_voxel;
tfce_stat_img.p_type = ['uncorrected TFCE p-values based on ' num2str(P) ' permutations'];

% Apply Storey's FDR
[~, TFCE_q_voxel_fdr, aprioriprob] = mafdr(TFCE_p_voxel);

% If aprioriprob > 0.99, fallback to BenjaminiHochberg
if aprioriprob > 0.99
    TFCE_p_voxel_fdr = mafdr(TFCE_p_voxel, 'BHFDR', true);
else
% Enforce constraint q >= p (as in SAS proc multtest)
    for j = 1:length(TFCE_q_voxel_fdr)
        if TFCE_q_voxel_fdr(j) < p_unc(j)
           TFCE_q_voxel_fdr(j) = p_unc(j);
        end
    end

end

tfce_stat_img_fdr = tfce_stat_img;
if exist('TFCE_p_voxel','var')
    tfce_stat_img_fdr.dat = TFCE_p_voxel_fdr;
    tfce_stat_img_fdr.dat_descrip = ['fdr corrected TFCE p-values based on ' num2str(P) ' permutations'];
    tfce_stat_img_fdr.p = TFCE_p_voxel_fdr;
    tfce_stat_img_fdr.p_type = ['fdr corrected TFCE p-values based on ' num2str(P) ' permutations'];
else
    tfce_stat_img_fdr.dat = TFCE_q_voxel_fdr;
    tfce_stat_img_fdr.dat_descrip = ['fdr corrected TFCE q-values based on ' num2str(P) ' permutations'];
    tfce_stat_img_fdr.p = TFCE_q_voxel_fdr;
    tfce_stat_img_fdr.p_type = ['fdr corrected TFCE q-values based on ' num2str(P) ' permutations'];
end

%% -------------------- Output --------------------
out.real_auc = real_auc;
out.real_t = real_t;
out.TFCE_real = TFCE_real;
out.TFCE_real_max = TFCE_real_max;
out.TFCE_null_max = maxTFCE_null;
out.p_TFCE_global = p_TFCE_global;
out.p_TFCE_voxel = TFCE_p_voxel;
out.tfce_stat_image = tfce_stat_img;
out.tfce_stat_image_fdr = tfce_stat_img_fdr;

fprintf('Global TFCE FWE p = %.4f\n',p_TFCE_global);
end

%% ========================================================================
% Helper: Classic TFCE (Smith & Nichols, 2009)
%% ========================================================================
function tfce = tfce_transform_3d(stat,H,E,conn)

if nargin<4, conn=26; end
tfce = zeros(size(stat),'single');

hvals = unique(stat(:));
hvals = hvals(hvals>0);
dh = diff([0; hvals]);

for i=1:numel(hvals)
   thr = hvals(i);
   bw = stat>=thr;
   CC = bwconncomp(bw,conn);
   for c=1:CC.NumObjects
       idx = CC.PixelIdxList{c};
       tfce(idx)=tfce(idx)+(thr^H)*(numel(idx)^E)*dh(i);
   end
end
end

%% ========================================================================
% Helper: AUC -> t (Hanley & McNeil, 1982)
%% ========================================================================
function t = auc_to_t(auc,n_pos,n_neg)

auc = max(min(auc,1-eps),eps);
Q1 = auc./(2-auc);
Q2 = (2*auc.^2)./(1+auc);

varAUC = (auc.*(1-auc) + (n_pos-1).*(Q1-auc.^2) + ...
         (n_neg-1).*(Q2-auc.^2))/(n_pos*n_neg);

SE = sqrt(max(varAUC,eps));
t = (auc-0.5)./SE;
end

function AUC = get_auc(dloc_template, dat_full, keep, alg, fixedC, Ccv)

   % --- dat ---
   dloc = dloc_template; % local copy
   dloc.dat = dat_full.dat(keep,:);

   % --- removed_voxels (LENGTH MUST EQUAL nvox) ---
   rv = true(dat_full.volInfo.nvox,1);
   rv(keep) = false;
   dloc.removed_voxels = rv;

   % --- volInfo fields ---
   % IMPORTANT: nvox must remain FULL MASK size
   dloc.volInfo.nvox = dat_full.volInfo.nvox;

   % In-mask voxel indices INTO FULL SPACE
   dloc.volInfo.wh_inmask = find(~rv);

   % xyzlist is ONLY in-mask coordinates
   dloc.volInfo.xyzlist = dat_full.volInfo.xyzlist(keep,:);

   % n_inmask must match dat rows
   dloc.volInfo.n_inmask = numel(keep);

   % --- predict ---
   if isempty(fixedC)
       if isempty(Ccv)
           [~,stats] = predict(dloc, ...
                               'algorithm_name',alg, ...
                               'verbose',0);
       else
           [~,stats] = predict(dloc, ...
                               'algorithm_name',alg, ...
                               'cv_assign',Ccv.Value, ...
                               'verbose',0);
       end
   else
       if isempty(Ccv)
           [~,stats] = predict(dloc, ...
                               'algorithm_name',alg, ...
                               'C',fixedC, ...
                               'verbose',0);
       else
           [~,stats] = predict(dloc, ...
                               'algorithm_name',alg, ...
                               'C',fixedC, ...
                               'cv_assign',Ccv.Value, ...
                               'verbose',0);
       end
   end

   % --- AUC ---
   scores = stats.dist_from_hyperplane_xval;
   ytrue = (dat_full.Y == 1);
   [~,~,~,AUC] = perfcurve(ytrue, scores, true);
end