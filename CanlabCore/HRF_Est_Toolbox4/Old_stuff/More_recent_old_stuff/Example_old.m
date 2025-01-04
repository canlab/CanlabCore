%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example code for estimating the HRF using the Inverse-Logit Model, a
% Finte Impluse Response Model and the Canonical HRF with 2 derivatives.
% Also the code illustrates our code for detecting model misspecification. 
%
% By Martin Lindquist and Tor Wager
% Edited  05/17/13
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load time course
%

mypath = which('ilogit');
if isempty(mypath), error('Cannot find directory with ilogit.m and other functions. Not on path?'); end
[mydir] = fileparts(mypath)

load(fullfile(mydir,'timecourse'))

tc = (tc- mean(tc))/std(tc);
len = length(tc);


%% Or: create your own
[xBF] = spm_get_bf(struct('dt', .5, 'name', 'hrf (with time and dispersion derivatives)', 'length', 32));
clear Xtrue
for i = 1:1, xx = conv(xBF.bf(:,i), [1 1 1 1 1 1 ]');
    Xtrue(:, i) = xx(1:66);
end
for i = 2:3, xx = conv(xBF.bf(:,i), [1 1]');
    Xtrue(:, i) = xx(1:66);
end
hrf = Xtrue * [1 .3 .2]';
xsecs = 0:.5:32;
hrf = hrf(1:length(xsecs));
hrf = hrf ./ max(hrf);
figure; plot(xsecs, hrf, 'k')
%hrf = hrf(1:4:end); % downsample to TR, if TR is > 0.5


R = randperm(640); R = sort(R(1:36));
Run = zeros(640,1);
for i=1:length(R), Run(R(i)) = 1; end;
true_sig = conv(Run, hrf);
true_sig = true_sig(1:640);

tc_noise = noise_arp(640, [.7 .2]);
tc = true_sig + 0 * tc_noise;
%figure; plot(tc);


%%

create_figure; subplot(3,1,1); han = plot(tc);
title('Sample time course'); drawnow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Settings
% 

TR = 0.5;
% T = round(30/TR);
% t = 1:T;                      % samples at which to get Logit HRF Estimate
T = 30;
FWHM = 4;                       % FWHM for residual scan
pval = 0.01;
df = 600;
alpha = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create stick function (sample event onsets)
% Variable R contains onset times
% Variable Run contains stick (a.k.a. delta or indicator) function

R = [3, 21, 56, 65, 109, 126, 163, 171, 216, 232, 269, 282, 323, 341, 376, 385, 429, 446, 483, 491, 536, 552, 589, 602];
Run = zeros(640,1);
for i=1:length(R), Run(R(i)) = 1; end;
Runc = {};
Runc{1} = Run;

try
    hold on;
    hh = plot_onsets(R,'k',-3,1);
    drawnow
catch
    disp('Couldn''t find function to add onset sticks to plot. Skipping.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using IL-function

% Choose mode (deterministic/stochastic)

mode = 1;   % 0 - deterministic aproach 
            % 1 - simulated annealing approach
            % Please note that when using simulated annealing approach you
            % may need to perform some tuning before use.

%[h1, fit1, e1, param] = Fit_Logit_allstim(tc,Run,t,mode);
V0 = [ 1 6 1 0.5 10 1 15];  
[VM, h1, fit1, e1, param] = Fit_Logit_allstim(V0,tc,TR, Runc, T);
[pv sres sres_ns1] = ResidScan(e1, FWHM);
[PowLoss1] = PowerLoss(e1, fit1, (len-7) , tc, TR, Runc, alpha);

hold on; han(2) = plot(fit1,'r');

disp('Summary: IL_function');

disp('Amplitude:'); disp(param(1));
disp('Time-to-peak:'); disp(param(2));
disp('Width:'); disp(param(3));

disp('MSE:'); disp((1/(len-1)*sum(e1.^2)));
disp('Mis-modeling:'); disp(pv);
disp('Power Loss:'); disp(PowLoss1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using FIR-model

% Choose mode (FIR/sFIR)

mode = 1;   % 0 - FIR 
            % 1 - smooth FIR
            
[h2, fit2, e2, param] = Fit_sFIR(tc,TR,Runc,T,mode);
[pv sres sres_ns2] = ResidScan(e2, FWHM);
[PowLoss2] = PowerLoss(e2, fit2, (len-T) , tc, TR, Runc, alpha);

hold on; han(3) = plot(fit2,'g');

disp('Summary: FIR');

disp('Amplitude'); disp(param(1));
disp('Time-to-peak'); disp(param(2));
disp('Width'); disp(param(3));

disp('MSE:'); disp((1/(len-1)*sum(e2.^2)));
disp('Mis-modeling'); disp(pv);
disp('Power Loss:'); disp(PowLoss2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using Canonical HRF + 2 derivatives

p=1;
      
[h3, fit3, e3, param, info] = Fit_Canonical_HRF(tc,TR,Runc,T,p);
[pv sres sres_ns3] = ResidScan(e3, FWHM);
[PowLoss3] = PowerLoss(e3, fit3, (len-p) , tc, TR, Runc, alpha);

hold on; han(4) = plot(fit3,'m');

legend(han,{'Data' 'IL' 'sFIR' 'DD'})


disp('Summary: Canonical + 2 derivatives');

disp('Amplitude'); disp(param(1));
disp('Time-to-peak'); disp(param(2));
disp('Width'); disp(param(3));

disp('MSE:'); disp((1/(len-1)*sum(e3.^2)));
disp('Mis-modeling'); disp(pv);
disp('Power Loss:'); disp(PowLoss3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure; 


subplot(3,2,5);
han2 = plot(h1,'r');
hold on; han2(2) = plot(h2,'g');
hold on; han2(3) = plot(h3,'m');
legend(han2,{'IL' 'sFIR' 'DD'})
title('Estimated HRF');


subplot(3,1,2); hold on;
hh = plot_onsets(R,'k',-3,1);
drawnow

han3 = plot(sres_ns1,'r');
hold on; han3(2) = plot(sres_ns2,'g');
hold on; han3(3) = plot(sres_ns3,'m');
hold on; plot((1:len),zeros(len,1),'--k');
legend(han3,{'IL' 'sFIR' 'DD'})
title('Mis-modeling (time course)');


subplot(3,2,6); hold on;

[s1] = Fit_sFIR(sres_ns1,TR,Runc,T,0);
[s2] = Fit_sFIR(sres_ns2,TR,Runc,T,0);
[s3] = Fit_sFIR(sres_ns3,TR,Runc,T,0);

han4 = plot(s1(1:T),'r');
hold on; han4(2) = plot(s2(1:T),'g');
hold on; han4(3) = plot(s3(1:T),'m');
hold on; plot((1:T),zeros(T,1),'--k');
legend(han4,{'IL' 'sFIR' 'DD'})
title('Mis-modeling (HRF)');
