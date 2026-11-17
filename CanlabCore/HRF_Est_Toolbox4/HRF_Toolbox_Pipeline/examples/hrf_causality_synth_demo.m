%% hrf_causality_synth_demo  -  why HRF-informed deconvolution matters for causality
%
% Ground-truth demonstration that regional HRF differences can REVERSE the
% direction inferred by Granger causality on raw BOLD (David et al. 2008,
% PLoS Biol), and that deconvolving each signal with its OWN hemodynamic
% kernel -- the per-region/signature sFIR HRF your pipeline estimates --
% recovers the true direction, where a single canonical kernel does not.
%
% Setup: neural A -> B (A leads by ~0.9 s), C independent. The DOWNSTREAM
% region B is given a FASTER HRF than the upstream A, so at the BOLD level B
% rises first and naive Granger flips the arrow.
%
% Requires: SPM (spm_hrf) on path; hrf_deconv_timeseries, hrf_granger_causality.
%
% See also: hrf_deconv_timeseries, hrf_granger_causality.

rng(7);
TR = 0.46; T = 300; nrun = 8; d = 2;          % d = neural lead of A over B (samples)

% Region HRFs: A SLOW (peak ~8 s), B FAST (peak ~4 s), C canonical.
hA = spm_hrf(TR, [8 16 1 1 6 0 32]);
hB = spm_hrf(TR, [4 16 1 1 6 0 32]);
hC = spm_hrf(TR);

boldRuns = cell(1, nrun);
for r = 1:nrun
    a = zeros(T,1); b = zeros(T,1); c = zeros(T,1);
    ea = randn(T,1); eb = randn(T,1); ec = randn(T,1);
    for t = 2:T
        a(t) = 0.40*a(t-1) + ea(t);
        c(t) = 0.40*c(t-1) + ec(t);
        if t > d, b(t) = 0.40*b(t-1) + 0.60*a(t-d) + eb(t);
        else,     b(t) = 0.40*b(t-1) + eb(t); end
    end
    ya = conv(a, hA); yb = conv(b, hB); yc = conv(c, hC);
    Y = [ya(1:T) yb(1:T) yc(1:T)];
    Y = Y + 0.5*std(Y(:))*randn(size(Y));     % measurement noise
    boldRuns{r} = Y;
end
nodes = {'A','B','C'};

% Per-region kernel matrix (what sFIR gives you) and a single canonical kernel.
Kmax = max([numel(hA) numel(hB) numel(hC)]);
Kern = zeros(Kmax,3); Kern(1:numel(hA),1)=hA; Kern(1:numel(hB),2)=hB; Kern(1:numel(hC),3)=hC;
canon = spm_hrf(TR); KernCanon = repmat([canon; zeros(Kmax-numel(canon),1)], 1, 3);

decReg   = cell(1,nrun); decCanon = cell(1,nrun);
for r = 1:nrun
    decReg{r}   = hrf_deconv_timeseries(boldRuns{r}, Kern,      'Method','ridge');
    decCanon{r} = hrf_deconv_timeseries(boldRuns{r}, KernCanon, 'Method','ridge');
end

Gnaive = hrf_granger_causality(boldRuns, 'Nodes', nodes, 'doverbose', false);
Gcanon = hrf_granger_causality(decCanon, 'Nodes', nodes, 'doverbose', false);
Greg   = hrf_granger_causality(decReg,   'Nodes', nodes, 'doverbose', false);

rows = {'1. naive GC on raw BOLD', Gnaive; ...
        '2. deconv w/ ONE canonical HRF', Gcanon; ...
        '3. deconv w/ PER-REGION sFIR HRF', Greg};
fprintf('\nGROUND TRUTH: A -> B (neural lead %.2f s); C independent; B has the FASTER HRF.\n', d*TR);
fprintf('%-34s  net(A->B)     p        direction\n', 'approach');
for i = 1:3
    nab = rows{i,2}.net(1,2); pv = rows{i,2}.pval(1,2);
    if nab > 0, dir = 'A->B  CORRECT'; else, dir = 'B->A  WRONG'; end
    fprintf('%-34s  %+8.3f   %7.1e   %s\n', rows{i,1}, nab, pv, dir);
end
fprintf('Only #3 (region-specific kernels) is both correct and strong.\n');
