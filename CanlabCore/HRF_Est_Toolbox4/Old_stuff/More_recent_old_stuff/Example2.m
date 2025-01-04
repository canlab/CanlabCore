%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example code for estimating the HRF using the Inverse-Logit Model, a
% Finte Impluse Response Model and the Canonical HRF with 2 derivatives.
% Also the code illustrates our code for detecting model misspecification. 
%
% This example illustrates the multi-stimulus version
%
% By Martin Lindquist
% Created: 05/26/10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load time course
%

mypath = which('ilogit');
if isempty(mypath), error('Cannot find directory with ilogit.m and other functions. Not on path?'); end
[mydir] = fileparts(mypath)


h = spm_hrf(1);
h = h./max(h);

RunA = zeros(60,10);
RunA(1,:) = 1;
RunA = reshape(RunA,600,1);
RunB = zeros(60,10);
RunB(31,:) = 1;
RunB = reshape(RunB,600,1);

tc = 2*conv(RunA,h) + conv(RunB,h);
tc = tc(1:600);
tc = tc+normrnd(0,1,600,1);

Run =[];
Run{1} = RunA;
Run{2} = RunB;

len = length(tc);

figure; subplot(3,1,1); han = plot(tc);
title('Sample time course'); drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Settings
% 

TR = 1;
T = round(30/TR);
t = 1:T;                        % samples at which to get Logit HRF Estimate
FWHM = 4;                       % FWHM for residual scan
pval = 0.01;
df = 600;
alpha = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create stick function (sample event onsets)
% Variable R contains onset times
% Variable Run contains stick (a.k.a. delta or indicator) function

RA = [1 61 121 181 241 301 361 421 481 541];
RB = [31 91 151 211 271 331 391 451 511 571];


try
    hold on;
    hh = plot_onsets(RA,'k',-3,1);
    hh = plot_onsets(RB,'r',-3,0.5);
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

     
            
            
[h1, fit1, e1, param] = Fit_Logit(tc,Run,t,mode);
[pv sres sres_ns1] = ResidScan(e1, FWHM);
[PowLoss1] = PowerLoss(e1, fit1, (len-7) , tc, TR, Run, alpha);

hold on; han(2) = plot(fit1,'r');

disp('Summary: IL_function');

disp('HRF - Event A');
disp('Amplitude:'); disp(param(1,1));
disp('Time-to-peak:'); disp(param(2,1));
disp('Width:'); disp(param(3,1));

disp('HRF - Event B');
disp('Amplitude:'); disp(param(1,2));
disp('Time-to-peak:'); disp(param(2,2));
disp('Width:'); disp(param(3,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using FIR-model

% Choose mode (FIR/sFIR)

mode = 1;   % 0 - FIR 
            % 1 - smooth FIR
            
[h2, fit2, e2, param] = Fit_sFIR(tc,TR,Run,T,mode);
[pv sres sres_ns2] = ResidScan(e2, FWHM);
[PowLoss2] = PowerLoss(e2, fit2, (len-T) , tc, TR, Run, alpha);

hold on; han(3) = plot(fit2,'g');

disp('Summary: FIR');

disp('HRF - Event A');
disp('Amplitude:'); disp(param(1,1));
disp('Time-to-peak:'); disp(param(2,1));
disp('Width:'); disp(param(3,1));

disp('HRF - Event B');
disp('Amplitude:'); disp(param(1,2));
disp('Time-to-peak:'); disp(param(2,2));
disp('Width:'); disp(param(3,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using Canonical HRF + 2 derivatives

p=1;
      
[h3, fit3, e3, param, info] = Fit_Canonical_HRF(tc,TR,Run,T,p);
[pv sres sres_ns3] = ResidScan(e3, FWHM);
[PowLoss3] = PowerLoss(e3, fit3, (len-p) , tc, TR, Run, alpha);

hold on; han(4) = plot(fit3,'m');

legend(han,{'Data' 'IL' 'sFIR' 'DD'})


disp('Summary: Canonical + 2 derivatives');

disp('HRF - Event A');
disp('Amplitude:'); disp(param(1,1));
disp('Time-to-peak:'); disp(param(2,1));
disp('Width:'); disp(param(3,1));

disp('HRF - Event B');
disp('Amplitude:'); disp(param(1,2));
disp('Time-to-peak:'); disp(param(2,2));
disp('Width:'); disp(param(3,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure; 


subplot(3,2,5);
han2 = plot(h1(:,1),'r');
hold on; han2(2) = plot(h2(:,1),'g');
hold on; han2(3) = plot(h3(:,1),'m');
legend(han2,{'IL' 'sFIR' 'DD'})
title('Estimated HRF - Event A');


subplot(3,1,2); hold on;
hh = plot_onsets(RA,'k',-3,1);
hh = plot_onsets(RB,'r',-3,0.5);
drawnow

han3 = plot(sres_ns1,'r');
hold on; han3(2) = plot(sres_ns2,'g');
hold on; han3(3) = plot(sres_ns3,'m');
hold on; plot((1:len),zeros(len,1),'--k');
legend(han3,{'IL' 'sFIR' 'DD'})
title('Mis-modeling (time course)');


subplot(3,2,6); hold on;

han4 = plot(h2(:,1),'r');
hold on; han4(2) = plot(h2(:,2),'g');
hold on; han4(3) = plot(h3(:,2),'m');
legend(han4,{'IL' 'sFIR' 'DD'})
title('Estimated HRF - Event A');
