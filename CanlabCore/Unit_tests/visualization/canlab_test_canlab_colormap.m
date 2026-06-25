function tests = canlab_test_canlab_colormap
%CANLAB_TEST_CANLAB_COLORMAP Unit tests for the central value->RGB mapping.
tests = functiontests(localfunctions);
end


% ---- solid -------------------------------------------------------------

function test_solid_is_constant(tc)
cm = canlab_colormap.solid([1 0 0]);
rgb = map(cm, [-3 0 2 9]');
tc.verifyEqual(rgb, repmat([1 0 0], 4, 1), 'AbsTol', 1e-12);
end


% ---- single ramp -------------------------------------------------------

function test_single_endpoints_and_midpoint(tc)
cm = canlab_colormap.single([1 0 0], [1 1 0], [0 4]);   % red -> yellow over 0..4
tc.verifyEqual(map(cm, 0),  [1 0 0],   'AbsTol', 1e-12);   % min colour at lo
tc.verifyEqual(map(cm, 4),  [1 1 0],   'AbsTol', 1e-12);   % max colour at hi
tc.verifyEqual(map(cm, 2),  [1 .5 0],  'AbsTol', 1e-12);   % midpoint
end

function test_single_clamps_outside_range(tc)
cm = canlab_colormap.single([1 0 0], [1 1 0], [0 4]);
tc.verifyEqual(map(cm, -10), [1 0 0], 'AbsTol', 1e-12);    % below -> min
tc.verifyEqual(map(cm,  99), [1 1 0], 'AbsTol', 1e-12);    % above -> max
end

function test_single_matches_render_blobs_formula(tc)
% render_blobs single map: w = clamp((v-lo)/(hi-lo)); rgb = w*maxcol + (1-w)*mincol.
mn = [1 0 0]; mx = [0 1 1]; rng = [-2 6];
cm = canlab_colormap.single(mn, mx, rng);
v = [-2 0 3 6]';
w = min(max((v - rng(1)) ./ (rng(2) - rng(1)), 0), 1);
expected = w .* mx + (1 - w) .* mn;
tc.verifyEqual(map(cm, v), expected, 'AbsTol', 1e-12);
end


% ---- split -------------------------------------------------------------

function test_split_positive_side(tc)
% positive: minpos (near 0) -> maxpos (extreme), over [posmin posmax].
cm = canlab_colormap.split([0 0 1], [0 1 1], [1 .5 0], [1 1 0], [-4 -2 2 4]);
tc.verifyEqual(map(cm, 2), [1 .5 0], 'AbsTol', 1e-12);   % minpos at posmin
tc.verifyEqual(map(cm, 4), [1 1 0],  'AbsTol', 1e-12);   % maxpos at posmax
tc.verifyEqual(map(cm, 3), [1 .75 0],'AbsTol', 1e-12);   % halfway
end

function test_split_negative_side(tc)
% negative: maxneg (near 0) -> minneg (extreme), over [negmax negmin].
cm = canlab_colormap.split([0 0 1], [0 1 1], [1 .5 0], [1 1 0], [-4 -2 2 4]);
tc.verifyEqual(map(cm, -2), [0 1 1], 'AbsTol', 1e-12);   % maxneg near zero
tc.verifyEqual(map(cm, -4), [0 0 1], 'AbsTol', 1e-12);   % minneg extreme
tc.verifyEqual(map(cm, -3), [0 .5 1],'AbsTol', 1e-12);   % halfway
end

function test_split_zero_is_uncoloured(tc)
cm = canlab_colormap.split([0 0 1], [0 1 1], [1 .5 0], [1 1 0], [-4 4]);
tc.verifyTrue(all(isnan(map(cm, 0))), 'value exactly 0 is uncoloured (NaN)');
end

function test_split_two_element_range_expands(tc)
% [lo hi] expands to [lo 0 0 hi]; positive uses 0..hi, negative lo..0.
cm = canlab_colormap.split([0 0 1], [0 1 1], [1 .5 0], [1 1 0], [-5 5]);
tc.verifyEqual(cm.range, [-5 0 0 5], 'AbsTol', 1e-12);
tc.verifyEqual(map(cm, 5),  [1 1 0], 'AbsTol', 1e-12);
tc.verifyEqual(map(cm, -5), [0 0 1], 'AbsTol', 1e-12);
end


% ---- indexed -----------------------------------------------------------

function test_indexed_rounds_and_clamps(tc)
cmap = [1 0 0; 0 1 0; 0 0 1];
cm = canlab_colormap.indexed(cmap);
tc.verifyEqual(map(cm, [1 2 3]'), cmap, 'AbsTol', 1e-12);
tc.verifyEqual(map(cm, 2.4), [0 1 0], 'AbsTol', 1e-12);   % rounds
tc.verifyEqual(map(cm, 99),  [0 0 1], 'AbsTol', 1e-12);   % clamps high
tc.verifyEqual(map(cm, -5),  [1 0 0], 'AbsTol', 1e-12);   % clamps low
end


% ---- from_render_args bridge ------------------------------------------

function test_from_render_args_types(tc)
s = canlab_colormap.from_render_args({'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}}, [-3 3]);
tc.verifyEqual(s.type, 'split');
tc.verifyEqual(s.range, [-3 0 0 3], 'AbsTol', 1e-12);

w = canlab_colormap.from_render_args({'maxcolor', [1 1 0], 'mincolor', [1 0 0]}, [0 5]);
tc.verifyEqual(w.type, 'single');
tc.verifyEqual(map(w, 5), [1 1 0], 'AbsTol', 1e-12);

so = canlab_colormap.from_render_args({'color', [1 0 0]}, []);
tc.verifyEqual(so.type, 'solid');

ix = canlab_colormap.from_render_args({'colormap', [1 0 0; 0 1 0]}, []);
tc.verifyEqual(ix.type, 'indexed');
end

function test_from_render_args_default_is_mango(tc)
% No colour spec -> the addblobs default (mango) split colormap.
d = canlab_colormap.from_render_args({}, [-2 2]);
tc.verifyEqual(d.type, 'split');
tc.verifyEqual(d.colors, {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]});
end


% ---- legend / lut ------------------------------------------------------

function test_legend_samples_split_spans_zero(tc)
% Split legend is ONE bar from extreme neg to extreme pos (through zero).
cm = canlab_colormap.split([0 0 1], [0 1 1], [1 .5 0], [1 1 0], [-4 4]);
[vals, rgb] = legend_samples(cm, 101);
tc.verifyEqual(vals(1),   -4, 'AbsTol', 1e-9, 'starts at extreme negative');
tc.verifyEqual(vals(end),  4, 'AbsTol', 1e-9, 'ends at extreme positive');
tc.verifyTrue(any(vals < 0) && any(vals > 0), 'spans zero');
tc.verifyEqual(rgb(1, :),   [0 0 1], 'AbsTol', 1e-9, 'extreme neg colour');
tc.verifyEqual(rgb(end, :), [1 1 0], 'AbsTol', 1e-9, 'extreme pos colour');
tc.verifyFalse(any(isnan(rgb(:))), 'no NaNs in legend (zero greyed)');
end

function test_lut_shape(tc)
cm = canlab_colormap.single([1 0 0], [1 1 0], [0 5]);
L = lut(cm, 256);
tc.verifyEqual(size(L), [256 3]);
tc.verifyEqual(L(1, :),   [1 0 0], 'AbsTol', 1e-9);
tc.verifyEqual(L(end, :), [1 1 0], 'AbsTol', 1e-9);
end


% ---- vectorized output -------------------------------------------------

function test_map_vector_shape(tc)
cm = canlab_colormap.split([0 0 1], [0 1 1], [1 .5 0], [1 1 0], [-4 4]);
rgb = map(cm, [-4 -2 2 4 0]');
tc.verifyEqual(size(rgb), [5 3]);
tc.verifyTrue(all(isnan(rgb(5, :))), 'zero row is NaN');
end


% ---- colorbar_ramp (continuous legend bar, no threshold gap) ------------

function test_colorbar_ramp_split_no_gap(tc)
% Split ramp: first half neg (minneg->maxneg), second half pos (minpos->maxpos),
% no greyed gap (unlike legend_samples). Ends are the extremes.
cm = canlab_colormap.split([0 0 1], [0 1 1], [1 .5 0], [1 1 0], [-4 4]);
r  = colorbar_ramp(cm, 64);
tc.verifyEqual(size(r), [64 3]);
tc.verifyFalse(any(isnan(r(:))), 'no gap/NaNs in colorbar ramp');
tc.verifyEqual(r(1, :),   [0 0 1], 'AbsTol', 1e-9, 'neg extreme = minneg');
tc.verifyEqual(r(end, :), [1 1 0], 'AbsTol', 1e-9, 'pos extreme = maxpos');
end

function test_colorbar_ramp_single_and_solid(tc)
cs = canlab_colormap.single([1 0 0], [1 1 0], [2 5]);
r  = colorbar_ramp(cs, 10);
tc.verifyEqual(r(1, :),   [1 0 0], 'AbsTol', 1e-9);
tc.verifyEqual(r(end, :), [1 1 0], 'AbsTol', 1e-9);
cm = canlab_colormap.solid([0.2 0.6 0.9]);
rs = colorbar_ramp(cm, 8);
tc.verifyEqual(rs, repmat([0.2 0.6 0.9], 8, 1), 'AbsTol', 1e-9, 'solid ramp is one colour');
end
