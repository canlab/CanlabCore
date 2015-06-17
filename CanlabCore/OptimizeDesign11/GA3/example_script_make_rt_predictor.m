cd /Users/tor/Documents/Tor_Documents/CurrentExperiments/Amygdala_fmri_class/Experiment_timing
load leftRhtSac
[X,d] = onsets2delta(sac,2);
figure;plot(X)
TR = 2;

% create fake RT data
rt = sac; for i = 1:2, rt{i} = 1000 + 100 * randn(size(sac{i}));,end

% build the model
[X,d,out] = rt2delta(sac,rt,TR);

% both build and plot
[X,d,out] = plotDesign(sac,rt,TR,2);


