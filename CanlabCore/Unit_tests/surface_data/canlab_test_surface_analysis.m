function tests = canlab_test_surface_analysis
% Unit tests for fmri_surface_data analysis methods (M6).
%
% Covers cat / horzcat, ttest, regress (native OLS), predict (delegated, weight
% map remapped to surface), and ica (delegated; skipped if the icatb toolbox is
% absent). Uses synthetic grayordinate objects -- no external toolbox required
% (except ica, which needs icatb_fastICA and is skipped otherwise).
%
% Run:    runtests('canlab_test_surface_analysis')
%
% :See also: cat, ttest, regress, predict, ica, fmri_surface_data

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
function setupOnce(t)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'Surface_tools'));
addpath(fullfile(here, '..', '..'));
assert(~isempty(which('fmri_surface_data')), 'fmri_surface_data not on path.');
t.TestData.figvis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
end

function teardownOnce(t)
close all force
set(0, 'DefaultFigureVisible', t.TestData.figvis);
end


% -------------------------------------------------------------------------
function test_cat_and_horzcat(t)
A = local_obj((1:20)' + (0:5)*0.3);     % 20 x 6
B = local_obj((1:20)' + (0:5)*0.3 + 1); % 20 x 6, same space

C = cat(A, B);
verifyEqual(t, size(C.dat), [20 12], 'cat should concatenate along maps.');
C2 = [A, B];
verifyEqual(t, C2.dat, C.dat, 'horzcat must match cat.');
verifyEqual(t, numel(C.removed_images), 12, 'removed_images must track maps.');

% Per-map fields concatenate
A.image_names = arrayfun(@(k) sprintf('a%d',k), 1:6, 'unif', 0)';
B.image_names = arrayfun(@(k) sprintf('b%d',k), 1:6, 'unif', 0)';
C = cat(A, B);
verifyEqual(t, numel(C.image_names), 12, 'image_names must concatenate.');

% Space mismatch errors
Adiff = A; Adiff.surface_space = 'fsaverage_164k';
verifyError(t, @() cat(A, Adiff), 'fmri_surface_data:cat:space');
end


% -------------------------------------------------------------------------
function test_ttest(t)
% Construct data with a known positive mean across maps.
base = ones(20, 1) * 2 + 0.01 * (1:20)';      % mostly positive
C = local_obj(base + 0.1 * randn(20, 10));
tobj = ttest(C);
verifyEqual(t, class(tobj), 'fmri_surface_data');
verifyEqual(t, size(tobj.dat), [20 1], 't-map is one column.');
verifyTrue(t, isfield(tobj.additional_info, 'statistic'), 'stats stored in additional_info.');
st = tobj.additional_info.statistic;
verifyEqual(t, numel(st.p), 20, 'p has one value per grayordinate.');
verifyTrue(t, all(tobj.dat > 0), 'Positive-mean data should give positive t.');
% Cross-check against MATLAB ttest at one grayordinate
[~, p1, ~, s1] = ttest(double(C.dat(5, :)));
% Note: fmri_data stores .dat as single, so the delegated t is single-precision.
verifyEqual(t, double(tobj.dat(5)), s1.tstat, 'RelTol', 1e-4, 't mismatch vs MATLAB ttest.');
verifyEqual(t, double(st.p(5)), p1, 'AbsTol', 1e-4, 'p mismatch vs MATLAB ttest.');
end


% -------------------------------------------------------------------------
function test_regress_ols(t)
n = 15;
X = [ones(n,1), (1:n)'/n];                    % intercept + linear
truebeta = [0.5; 3];
C = local_obj((1:20)'*0 + (X * truebeta)' + 0.001*randn(20, n));  % 20 x n, ~linear in X
b = regress(C, X);
verifyEqual(t, size(b.dat), [20 2], 'betas: one column per regressor.');
% Cross-check coefficients against MATLAB backslash at one grayordinate
g = 9;
bm = X \ double(C.dat(g, :))';
verifyEqual(t, double(b.dat(g, :))', bm, 'AbsTol', 1e-5, 'OLS betas mismatch.');
st = b.additional_info.statistic;
verifyEqual(t, size(st.t), [20 2], 't-stats per beta.');
verifyEqual(t, st.dfe, n - 2, 'wrong residual df.');
% Recovered slope should be near the true slope (.dat is single)
verifyEqual(t, double(mean(b.dat(:,2))), truebeta(2), 'AbsTol', 0.05, 'slope not recovered.');
% Errors with no design
verifyError(t, @() regress(local_obj(randn(20,5))), 'fmri_surface_data:regress:noX');
end


% -------------------------------------------------------------------------
function test_predict_delegation(t)
% Outcome linearly related to a latent feature -> cv prediction should work and
% return a surface weight map of the right size.
n = 30;
w = randn(20, 1);
Y = (1:n)';
D = w * (Y' - mean(Y)) + randn(20, n);        % each grayordinate ~ scaled Y
C = local_obj(D);
C.Y = Y;
[cverr, stats] = predict(C, 'algorithm_name', 'cv_lassopcr', 'nfolds', 3, 'verbose', 0);
verifyTrue(t, isfinite(cverr), 'cverr should be finite.');
verifyEqual(t, class(stats.weight_obj), 'fmri_surface_data', ...
    'weight_obj must be remapped to a surface object.');
verifyEqual(t, size(stats.weight_obj.dat, 1), 20, 'weight map has one row per grayordinate.');
verifyEqual(t, numel(stats.yfit), n, 'one cross-validated fit per observation.');
% No-Y guard
Cn = local_obj(randn(20, n));
verifyError(t, @() predict(Cn), 'fmri_surface_data:predict:noY');
end


% -------------------------------------------------------------------------
function test_ica_if_available(t)
if isempty(which('icatb_fastICA'))
    t.assumeFail('icatb_fastICA (GIFT/icatb toolbox) not on path; ica delegation skipped.');
end
n = 500; nm = 30;
S = randn(3, nm); A = randn(n, 3); D = A * S + 0.1 * randn(n, nm);
o = local_obj_n(D);
ic = ica(o, 3);
verifyEqual(t, class(ic), 'fmri_surface_data');
verifyEqual(t, size(ic.dat, 1), n, 'IC maps have one row per grayordinate.');
end


% =========================================================================
function o = local_obj(D)
% 20-grayordinate cortex-only object (10 verts per hemisphere).
mL = struct('struct','CORTEX_LEFT','type','surf','start',1,'count',10,'numvert',32492,'vertlist',0:9,'voxlist',[]);
mR = struct('struct','CORTEX_RIGHT','type','surf','start',11,'count',10,'numvert',32492,'vertlist',0:9,'voxlist',[]);
bm = struct('type','dense','length',20,'models',{{mL,mR}},'vol',[]);
bm.grayordinate_type = 'cortex_only'; bm.cluster = [];
o = fmri_surface_data('dat', single(D), 'brain_model', bm, ...
    'surface_space', 'fsLR_32k', 'intent', 'dscalar');
end

function o = local_obj_n(D)
% n-grayordinate cortex-only object (n/2 verts per hemisphere).
n = size(D,1); h = n/2;
mL = struct('struct','CORTEX_LEFT','type','surf','start',1,'count',h,'numvert',32492,'vertlist',0:h-1,'voxlist',[]);
mR = struct('struct','CORTEX_RIGHT','type','surf','start',h+1,'count',h,'numvert',32492,'vertlist',0:h-1,'voxlist',[]);
bm = struct('type','dense','length',n,'models',{{mL,mR}},'vol',[]);
bm.grayordinate_type = 'cortex_only'; bm.cluster = [];
o = fmri_surface_data('dat', single(D), 'brain_model', bm, ...
    'surface_space', 'fsLR_32k', 'intent', 'dscalar');
end
