%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example code for estimating the HRF using the Inverse-Logit Model, a
% Finte Impluse Response Model and the Canonical HRF with 2 derivatives.
% Also the code illustrates our code for detecting model misspecification. 
%
% By Martin Lindquist and Tor Wager
% Created  03/09/23
% Last edited 03/09/23
% Added multi-simulus simulation to Example.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create time course
%

mypath = which('ilogit');
if isempty(mypath), error('Cannot find directory with ilogit.m and other functions. Not on path?'); end
[mydir] = fileparts(mypath);

load(fullfile(mydir,'timecourse'))

tc = (tc- mean(tc))/std(tc);
len = length(tc);


%%
% Create HRFs

[xBF] = spm_get_bf(struct('dt', .5, 'name', 'hrf (with time and dispersion derivatives)', 'length', 32));
clear Xtrue1
for i = 1:1, xx = conv(xBF.bf(:,i), [1 1 1 1 1 1 ]');
    Xtrue1(:, i) = xx(1:66);
end
for i = 2:3, xx = conv(xBF.bf(:,i), [1 1]');
    Xtrue1(:, i) = xx(1:66);
end
hrf1 = Xtrue1 * [1 .5 .3]';


clear Xtrue2
for i = 1:1, xx = conv(xBF.bf(:,i), [1 1 ]');
    Xtrue2(:, i) = xx(1:66);
end
for i = 2:3, xx = conv(xBF.bf(:,i), [1 1]');
    Xtrue2(:, i) = xx(1:66);
end
hrf2 = Xtrue2 * [1 0 0]';

xsecs = 0:.5:32;

hrf1 = [ 0; 0; hrf1];
hrf1 = hrf1(1:length(xsecs));
hrf1 = hrf1 ./ max(hrf1);

hrf2 = [ 0; 0; hrf2];
hrf2 = hrf2(1:length(xsecs));
hrf2 = hrf2 ./ max(hrf2);

figure; plot(xsecs, hrf1, 'k')
hold; plot(xsecs, hrf2, 'g')


%hrf = hrf(1:4:end); % downsample to TR, if TR is > 0.5

% Create stimuli

R = randperm(640); 
t1 = 1:18;
t2 = 19:36;
R1 = sort(R(t1));
R2 = sort(R(t2));

Run1 = zeros(640,1);
Run2 = zeros(640,1);
for i=1:length(R1), Run1(R1(i)) = 1; Run2(R2(i)) = 1; end;

% Create timecourse

beta1 = 1; beta2 = 0.8;
true_sig = beta1*conv(Run1, hrf1) + beta2*conv(Run2, hrf2);
true_sig = true_sig(1:640);

tc_noise = noise_arp(640, [.3 0]);
tc = true_sig +  0.5 * tc_noise;
% tc = true_sig;
%figure; plot(tc);


Runc{1} = Run1;
Runc{2} = Run2;


%%

create_figure; subplot(3,1,1); han = plot(tc);
title('Sample time course'); drawnow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Settings
% 

TR = 0.5;
%T = round(30/TR);
T = 30;
t = 1:TR:T;                        % samples at which to get Logit HRF Estimate
FWHM = 4;                       % FWHM for residual scan
pval = 0.01;
df = 600;
alpha = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create stick function (sample event onsets)
% Variable R contains onset times
% Variable Run contains stick (a.k.a. delta or indicator) function

% R = [3, 21, 56, 65, 109, 126, 163, 171, 216, 232, 269, 282, 323, 341, 376, 385, 429, 446, 483, 491, 536, 552, 589, 602];
% Run = zeros(640,1);
% for i=1:length(R), Run(R(i)) = 1; end;
% 

try
    hold on;
    hh = plot_onsets(R1,'r',-3,1, 1);
    hh = plot_onsets(R2,'g',-3,1, 1);
    drawnow
catch
    disp('Couldn''t find function to add onset sticks to plot. Skipping.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using IL-function

% Choose mode (deterministic/stochastic)

mode = 0;   % 0 - deterministic aproach 
            % 1 - simulated annealing approach
            % Please note that when using simulated annealing approach you
            % may need to perform some tuning before use.

[h1, fit1, e1, param] = Fit_Logit2(tc,TR,Runc,T,mode);
[pv sres sres_ns1] = ResidScan(e1, FWHM);
[PowLoss1] = PowerLoss(e1, fit1, (len-7) , tc, TR, Runc, alpha);

hold on; han(2) = plot(fit1,'r');

disp('Summary: IL_function');

disp('Amplitude:'); disp(param(1,:));
disp('Time-to-peak:'); disp(param(2,:)*TR);
disp('Width:'); disp(param(3,:)*TR);

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

disp('Amplitude'); disp(param(1,:));
disp('Time-to-peak'); disp(param(2,:)*TR);
disp('Width'); disp(param(3,:)*TR);

disp('MSE:'); disp((1/(len-1)*sum(e2.^2)));
disp('Mis-modeling'); disp(pv);
disp('Power Loss:'); disp(PowLoss2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using Canonical HRF + 2 derivatives

p=1;
      
[h3, fit3, e3, param, info] = Fit_Canonical_HRF(tc,TR,Runc,30,p);
[pv sres sres_ns3] = ResidScan(e3, FWHM);
[PowLoss3] = PowerLoss(e3, fit3, (len-p) , tc, TR, Runc, alpha);

hold on; han(4) = plot(fit3,'m');


disp('Summary: Canonical + 2 derivatives');

disp('Amplitude'); disp(param(1,:));
disp('Time-to-peak'); disp(param(2,:)*TR);
disp('Width'); disp(param(3,:)*TR);

disp('MSE:'); disp((1/(len-1)*sum(e3.^2)));
disp('Mis-modeling'); disp(pv);
disp('Power Loss:'); disp(PowLoss3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using Bsplines

[h4, fit4, e4, param] = Fit_Spline(tc, TR, Runc, T);

[pv sres sres_ns4] = ResidScan(e4, FWHM);
[PowLoss4] = PowerLoss(e4, fit4, (len-p) , tc, TR, Runc, alpha);

hold on; han(5) = plot(fit4,'b');


disp('Summary: Spline');

disp('Amplitude'); disp(param(1,:));
disp('Time-to-peak'); disp(param(2,:)*TR);
disp('Width'); disp(param(3,:)*TR);

disp('MSE:'); disp((1/(len-1)*sum(e4.^2)));
disp('Mis-modeling'); disp(pv);
disp('Power Loss:'); disp(PowLoss4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using non-linear gamma

[h5, fit5, e5, param] = Fit_NLgamma(tc, TR, Runc, T);

[pv sres sres_ns5] = ResidScan(e5, FWHM);
[PowLoss5] = PowerLoss(e5, fit5, (len-p) , tc, TR, Runc, alpha);

hold on; han(6) = plot(fit5,'b');

legend(han,{'Data' 'IL' 'sFIR' 'DD' 'Spline' 'NL'})


disp('Summary: NL gamma');

disp('Amplitude'); disp(param(1,:));
disp('Time-to-peak'); disp(param(2,:)*TR);
disp('Width'); disp(param(3,:)*TR);

disp('MSE:'); disp((1/(len-1)*sum(e5.^2)));
disp('Mis-modeling'); disp(pv);
disp('Power Loss:'); disp(PowLoss5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure; 
%%

subplot(3,2,5); hold on;
plot(xsecs, beta1*hrf1, 'k')
plot(xsecs, beta2*hrf2, 'k')

xsecs1 = xsecs(1:length(h1));
han2 = plot(xsecs1, h1,'r');

xsecs2 = xsecs(1:length(h2));
han2(3:4) = plot(xsecs2, h2,'g');
xsecs3 = xsecs(1:length(h3));
han2(5:6) = plot(xsecs3, h3,'m');
xsecs4 = xsecs(1:length(h4));
han2(7:8) = plot(xsecs4, h4,'b');
xsecs5 = xsecs(1:length(h5));
han2(9:10) = plot(xsecs5, h5,'y');
legend(han2,{'IL1' 'IL2' 'sFIR1' 'sFIR2' 'DD1' 'DD2' 'Spline1' 'Spline2' 'NL1' 'NL2'})

title('Estimated HRF');


subplot(3,1,2); hold on;
hh = plot_onsets(R1,'r',-3,1, 1);
hh = plot_onsets(R2,'g',-3,1, 1);
drawnow

han3 = plot(sres_ns1,'r');
hold on; han3(2) = plot(sres_ns2,'g');
hold on; han3(3) = plot(sres_ns3,'m');
hold on; han3(4) = plot(sres_ns4,'b');
hold on; han3(5) = plot(sres_ns5,'y');

hold on; plot((1:len),zeros(len,1),'--k');
legend(han3,{'IL' 'sFIR' 'DD' 'Spline' 'NL'})
title('Mis-modeling (time course)');


subplot(3,2,6); hold on;

[s1] = Fit_sFIR(sres_ns1,TR,Runc,T,0);
[s2] = Fit_sFIR(sres_ns2,TR,Runc,T,0);
[s3] = Fit_sFIR(sres_ns3,TR,Runc,T,0);
[s4] = Fit_sFIR(sres_ns4,TR,Runc,T,0);
[s5] = Fit_sFIR(sres_ns5,TR,Runc,T,0);

han4 = plot(s1(1:T),'r');
hold on; han4(2) = plot(s2(1:T),'g');
hold on; han4(3) = plot(s3(1:T),'m');
hold on; han4(4) = plot(s4(1:T),'b');
hold on; han4(5) = plot(s5(1:T),'y');
hold on; plot((1:T),zeros(T,1),'--k');

legend(han4,{'IL' 'sFIR' 'DD' 'Spline' 'NL'})
title('Mis-modeling (HRF)');

%%


figure; hold on;

% plot(xsecs, beta1*hrf1, 'k--', 'LineWidth', 2)
% plot(xsecs, beta2*hrf2, 'k--', 'LineWidth', 2)
% 
% xsecs1 = xsecs(1:length(h1));
% han2 = plot(xsecs1, h1,'r', 'LineWidth', 2);
% xsecs2 = xsecs(1:length(h2));
% han2(3:4) = plot(xsecs2, h2,'g', 'LineWidth', 2);
% xsecs3 = xsecs(1:length(h3));
% han2(5:6) = plot(xsecs3, h3,'m', 'LineWidth', 2);
% xsecs4 = xsecs(1:length(h4));
% han2(7:8) = plot(xsecs4, h4,'b', 'LineWidth', 2);
% xsecs5 = xsecs(1:length(h5));
% han2(9:10) = plot(xsecs5, h5,'y', 'LineWidth', 2);
% 
% legend(han2,{'IL1' 'IL2' 'sFIR1' 'sFIR2' 'DD1' 'DD2' 'Spline1' 'Spline2' 'NL1' 'NL2'})
% title('Estimated HRF');

subplot 121
hold
plot(xsecs, beta1*hrf1, 'k--', 'LineWidth', 2)

xsecs1 = xsecs(1:length(h1));
han2 = plot(xsecs1, h1(:,1),'r', 'LineWidth', 2);
xsecs2 = xsecs(1:length(h2));
han2(2) = plot(xsecs2, h2(:,1),'g', 'LineWidth', 2);
xsecs3 = xsecs(1:length(h3));
han2(3) = plot(xsecs3, h3(:,1),'m', 'LineWidth', 2);
xsecs4 = xsecs(1:length(h4));
han2(4) = plot(xsecs4, h4(:,1),'b', 'LineWidth', 2);
xsecs5 = xsecs(1:length(h5));
han2(5) = plot(xsecs5, h5(:,1),'y', 'LineWidth', 2);

legend(han2,{'IL1' 'sFIR1' 'DD1' 'Spline1' 'NL1'})
title('Estimated HRF');

subplot 122
hold
plot(xsecs, beta2*hrf2, 'k--', 'LineWidth', 2)

xsecs1 = xsecs(1:length(h1));
han2 = plot(xsecs1, h1(:,2),'r', 'LineWidth', 2);
xsecs2 = xsecs(1:length(h2));
han2(2) = plot(xsecs2, h2(:,2),'g', 'LineWidth', 2);
xsecs3 = xsecs(1:length(h3));
han2(3) = plot(xsecs3, h3(:,2),'m', 'LineWidth', 2);
xsecs4 = xsecs(1:length(h4));
han2(4) = plot(xsecs4, h4(:,2),'b', 'LineWidth', 2);
xsecs5 = xsecs(1:length(h5));
han2(5) = plot(xsecs5, h5(:,2),'y', 'LineWidth', 2);

legend(han2,{'IL2' 'sFIR2' 'DD2' 'Spline2' 'NL2'})
title('Estimated HRF');

