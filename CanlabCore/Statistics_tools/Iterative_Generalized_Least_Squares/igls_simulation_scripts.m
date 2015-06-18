%%  LITTLE dummy setup test
    len = 200; sub = 20;
    x = zeros(len,sub);
    x(11:20,:) = 2;                   % create signal
    x(111:120,:) = 2;                   % create signal
    c = normrnd(0.5,1,sub,1);       % slope between-subjects variations
    d = normrnd(0.5,1,sub,1);         % intercept between-subjects variations
    % Add between-subjects error (random effects) and measurement noise
    % (within-subjects error)
    sig = 0.1;
    for i=1:sub
        y(:,i) = d(i) + c(i).*x(:,i) + normrnd(0,sig,len,1);
    end
    out = igls(y, x, 'type', 'i');  % for igls
    disp('Input random-effect variances: '); disp(std([d c]))
    disp('Est.  random-effect variances: '); disp(sqrt(out.betastar)');
  
%%  VERY LITTLE dummy setup test
    len = 10; sub = 5;
    x = zeros(len,sub);
    x(1:2:10,:) = 2;                   % create signal
    c = normrnd(0.5,1,sub,1);       % slope between-subjects variations
    d = normrnd(0.5,1,sub,1);         % intercept between-subjects variations
    % Add between-subjects error (random effects) and measurement noise
    % (within-subjects error)
    clear y
    sig = 1;
    for i=1:sub
        y(:,i) = d(i) + c(i).*x(:,i) + noise_arp(len, [.5]);  %normrnd(0,sig,len,1);
    end
    out = igls(y, x, 'type', 'i', 'ar', 1);  % for igls
    disp('Input random-effect variances: '); disp(std([d c]))
    disp('Est.  random-effect variances: '); disp(sqrt(out.betastar)');
    
    
%% Basic test of igls variance param p-values for getting false pos. rates and power
clear *pvals* *fpr*

% SETUP SIM. INPUT PARAMS
mytype = 'i';
iterations = 1000;


len = 100; sub = 10;

fixpvals = zeros(iterations, 2);
randpvals = zeros(iterations, 2);
LRTrandpvals = zeros(iterations, 2);
chic = zeros(iterations, 1);
chid = zeros(iterations, 1);

x = zeros(len,sub);
y = x;
x(1:2:end,:) = 1;                   % create signal

noise_std = 2.0;  % within-subjects noise std

% ------------------------------------
%
% False pos rate
%
% ------------------------------------

pvals = zeros(iterations, 2);

c_fixed = 0;  % slope pop. average
d_fixed = 0;   % intercept pop. average
c_rand = 0;    % slope std across Ss
d_rand = 0;    % intercept std.

SUM = [];
SUM.input_options = struct('iterations', iterations, 'mytype', mytype, 'len', len, 'sub', sub, 'noise_std', noise_std, 'x', x);

% --------------------------------------------------------------------------
% This block: No true random effects, so test false positive rate for
% p_randvariance

verbstr = 'verbose';

for i = 1:iterations
  
    c = normrnd(c_fixed,c_rand,sub,1);       % slope between-subjects variations
    d = normrnd(d_fixed,d_rand,sub,1);         % intercept between-subjects variations
    
    % Add between-subjects error (random effects) and measurement noise
    % (within-subjects error)
    y = zeros(len, sub);
    for s=1:sub
        y(:,s) = d(s) + c(s) .* x(:,s) + normrnd(0,noise_std,len,1);
    end
    out = igls(y, x, 'type', mytype, verbstr);  % for igls
 
    fpr_converged(i) = out.isconverged;
    fpr_betastar(i, :) = out.betastar';
    
    pvals(i, :) = out.p_randvariance';

    fixpvals(i, :) = out.p;
    randpvals(i, :) = out.p_randvariance'; %[out.p_randvariance_d out.p_randvariance_c];
    LRTrandpvals(i, :) = out.pLRT_randvariance';
    
    
    if mod(i, 10) == 0, fprintf(1, '%3.0f ', i); end
    
    verbstr = 'noverbose';
end

pvals_fpr = pvals;
SUM.fixpvals_fpr = fixpvals;
SUM.randpvals_fpr = randpvals;
SUM.LRTrandpvals_fpr = LRTrandpvals;

SUM.fixfpr =  sum(SUM.fixpvals_fpr < .05) ./ iterations;
SUM.randfpr =  sum(SUM.LRTrandpvals_fpr < .05) ./ iterations;
fprintf('Fixed FPR: %3.5f\t%3.5f\tRFX FPR: %3.5f\t%3.5f\n', SUM.fixfpr, SUM.randfpr);

% ------------------------------------
%
% Power
%
% ------------------------------------

pvals = zeros(iterations, 2);
fixpvals = zeros(iterations, 2);
randpvals = zeros(iterations, 2);
LRTrandpvals = zeros(iterations, 2);

c_fixed = .5;  % slope pop. average
d_fixed = .5;   % intercept pop. average
c_rand = 0.5;    % slope std across Ss
d_rand = 0.5;    % intercept std.
% % % noise_std = 1;  % within-subjects noise std

SUM.input_options.c_fixed = c_fixed;
SUM.input_options.d_fixed = d_fixed;
SUM.input_options.c_rand = c_rand;
SUM.input_options.d_rand = d_rand;

verbstr = 'verbose';

for i = 1:iterations
  
    c = normrnd(c_fixed,c_rand,sub,1);       % slope between-subjects variations
    d = normrnd(d_fixed,d_rand,sub,1);         % intercept between-subjects variations
    
    % Add between-subjects error (random effects) and measurement noise
    % (within-subjects error)
    for s=1:sub
        y(:,s) = d(s) + c(s) .* x(:,s) + normrnd(0,noise_std,len,1);
    end
    out = igls(y, x, 'type', mytype, 'noverbose');  % for igls
 
    pow_converged(i) = out.isconverged;
    pow_betastar(i, :) = out.betastar';
    pvals(i, :) = out.p_randvariance';
    
    fixpvals(i, :) = out.p;
    randpvals(i, :) = out.p_randvariance'; %[out.p_randvariance_d out.p_randvariance_c];
    LRTrandpvals(i, :) = out.pLRT_randvariance';
    
    if mod(i, 10) == 0, fprintf(1, '%3.0f ', i); end
    
        verbstr = 'noverbose';
end

SUM.pvals_pow = pvals;
SUM.fixpvals_pow = fixpvals;
SUM.randpvals_pow = randpvals;
SUM.LRTrandpvals_pow = LRTrandpvals;

pow =  sum(pvals < .05) ./ iterations;
fprintf('Power: %3.5f\t%3.5f\n', pow);


% ------------------------------------
%
% ROC plot
%
% ------------------------------------

create_figure('IGLS ROC plot');
clear fpr pow fix* rand* LRT*
thresholds = .001:.001:.5;
for pthr = 1:length(thresholds)
%     fpr(pthr) =  sum(SUM.pvals_fpr(:,2) < thresholds(pthr)) ./ iterations;
%     pow(pthr) =  sum(SUM.pvals_pow(:,2) < thresholds(pthr)) ./ iterations;
    
    fix_fpr(pthr) =  sum(SUM.fixpvals_fpr(:,2) < thresholds(pthr)) ./ iterations;
    fix_pow(pthr) =  sum(SUM.fixpvals_pow(:,2) < thresholds(pthr)) ./ iterations;
    
    rand_fpr(pthr) =  sum(SUM.randpvals_fpr(:,2) < thresholds(pthr)) ./ iterations;
    rand_pow(pthr) =  sum(SUM.randpvals_pow(:,2) < thresholds(pthr)) ./ iterations;
    
    LRTrand_fpr(pthr) =  sum(SUM.LRTrandpvals_fpr(:,2) < thresholds(pthr)) ./ iterations;
    LRTrand_pow(pthr) =  sum(SUM.LRTrandpvals_pow(:,2) < thresholds(pthr)) ./ iterations;
       
end
 hold on;
 clear h

 scatter(fix_fpr, fix_pow, 34, log(thresholds), 'filled');
h(1) = plot(sort(fix_fpr), sort(fix_pow), 'k');

scatter(LRTrand_fpr, LRTrand_pow, 34, log(thresholds), 'filled');
h(2) = plot(sort(LRTrand_fpr), sort(LRTrand_pow), 'k:');

wh = find(thresholds - .05 == min(abs(thresholds - .05)));
plot([fix_fpr(wh) fix_fpr(wh)], [0 fix_pow(wh)], 'k--');
text(fix_fpr(wh) + .02, .1, 'p < .05','FontSize', 18);

wh = find(thresholds - .001 == min(abs(thresholds - .001)));
plot([fix_fpr(wh) fix_fpr(wh)], [0 fix_pow(wh)], 'k--');
text(fix_fpr(wh) + .01, .1, 'p < .001','FontSize', 18);

ylabel('Power'); xlabel('Actual false positive rate')
legend(h, {'Fixed effects' 'Random effects'});
set(gca, 'XLim', [-.01 .1]);



%% TEST COVARIATE


%% Basic test of igls variance param p-values for getting false pos. rates and power
clear *pvals* *fpr*

% SETUP SIM. INPUT PARAMS
mytype = 'i';
iterations = 1000;


len = 100; sub = 10;

fixpvals = zeros(iterations, 3);
randpvals = zeros(iterations, 2);
LRTrandpvals = zeros(iterations, 2);

x = zeros(len,sub);
y = x;
x(1:2:end,:) = 1;                   % create signal

noise_std = 2.0;  % within-subjects noise std

% ------------------------------------
%
% False pos rate
%
% ------------------------------------

pvals = zeros(iterations, 3);

c_fixed = 0;  % slope pop. average
d_fixed = 0;   % intercept pop. average
c_rand = 0;    % slope std across Ss
d_rand = 0;    % intercept std.
c_cov = 0;      % slope x 2nd-level covariate

SUM = [];
SUM.input_options = struct('iterations', iterations, 'mytype', mytype, 'len', len, 'sub', sub, 'noise_std', noise_std, 'x', x, 'c_cov', c_cov);

% --------------------------------------------------------------------------
% This block: No true random effects, so test false positive rate for
% p_randvariance

verbstr = 'verbose';

for i = 1:iterations
  
    cov_vals = randn(sub, 1); %+ sqrt(c_cov) .* randn(sub, 1);
    
    c = c_cov .* cov_vals + normrnd(c_fixed,c_rand,sub,1);       % slope between-subjects variations
    d = normrnd(d_fixed,d_rand,sub,1);                          % intercept between-subjects variations
    
    % Add between-subjects error (random effects) and measurement noise
    % (within-subjects error)
    y = zeros(len, sub);
    for s=1:sub
        y(:,s) = d(s) + c(s) .* x(:,s) + normrnd(0,noise_std,len,1);
    end
    out = igls(y, x, 'type', mytype, verbstr, 'covariate', cov_vals);  % for igls
 
    fpr_converged(i) = out.isconverged;
    fpr_betastar(i, :) = out.betastar';
    
    fixpvals(i, :) = out.p;
    randpvals(i, :) = out.p_randvariance'; %[out.p_randvariance_d out.p_randvariance_c];
    LRTrandpvals(i, :) = out.pLRT_randvariance';
    
    
    if mod(i, 10) == 0, fprintf(1, '%3.0f ', i); end
    
    verbstr = 'noverbose';
end

SUM.fixpvals_fpr = fixpvals;
SUM.randpvals_fpr = randpvals;
SUM.LRTrandpvals_fpr = LRTrandpvals;

SUM.fixfpr =  sum(SUM.fixpvals_fpr < .05) ./ iterations;
SUM.randfpr =  sum(SUM.LRTrandpvals_fpr < .05) ./ iterations;
fprintf('Fixed FPR: %3.5f\t%3.5f\tRFX FPR: %3.5f\t%3.5f\n', SUM.fixfpr, SUM.randfpr);

% ------------------------------------
%
% Power
%
% ------------------------------------

fixpvals = zeros(iterations, 3);
randpvals = zeros(iterations, 2);
LRTrandpvals = zeros(iterations, 2);

c_fixed = .5;  % slope pop. average
d_fixed = .5;   % intercept pop. average
c_rand = 0.15;    % slope std across Ss
d_rand = 0.15;    % intercept std.
c_cov = 0.25;      % slope x 2nd-level covariate

SUM.input_options.c_fixed = c_fixed;
SUM.input_options.d_fixed = d_fixed;
SUM.input_options.c_rand = c_rand;
SUM.input_options.d_rand = d_rand;
SUM.input_options.c_cov = c_cov;

verbstr = 'verbose';

for i = 1:iterations
  
    cov_vals = randn(sub, 1); %+ sqrt(c_cov) .* randn(sub, 1);
    
    c = c_cov .* cov_vals + normrnd(c_fixed,c_rand,sub,1);       % slope between-subjects variations
    d = normrnd(d_fixed,d_rand,sub,1);         % intercept between-subjects variations
    
    % Add between-subjects error (random effects) and measurement noise
    % (within-subjects error)
    y = zeros(len, sub);
    for s=1:sub
        y(:,s) = d(s) + c(s) .* x(:,s) + normrnd(0,noise_std,len,1);
    end
    
    out = igls(y, x, 'type', mytype, verbstr, 'covariate', cov_vals);  % for igls
 
    pow_converged(i) = out.isconverged;
    pow_betastar(i, :) = out.betastar';
    
    fixpvals(i, :) = out.p;
    randpvals(i, :) = out.p_randvariance'; %[out.p_randvariance_d out.p_randvariance_c];
    LRTrandpvals(i, :) = out.pLRT_randvariance';
    
    if mod(i, 10) == 0, fprintf(1, '%3.0f ', i); end
    
        verbstr = 'noverbose';
end

SUM.fixpvals_pow = fixpvals;
SUM.randpvals_pow = randpvals;
SUM.LRTrandpvals_pow = LRTrandpvals;

SUM.fixpow =  sum(SUM.fixpvals_pow < .05) ./ iterations;
SUM.randpow =  sum(SUM.LRTrandpvals_pow < .05) ./ iterations;
fprintf('Fixed FPR: %3.5f\t%3.5f\tRFX FPR: %3.5f\t%3.5f\n', SUM.fixpow, SUM.randpow);


% ------------------------------------
%
% ROC plot
%
% ------------------------------------
%
create_figure('IGLS ROC plot');
clear fpr pow fix* rand* LRT*
thresholds = .001:.001:.5;
for pthr = 1:length(thresholds)
%     fpr(pthr) =  sum(SUM.pvals_fpr(:,2) < thresholds(pthr)) ./ iterations;
%     pow(pthr) =  sum(SUM.pvals_pow(:,2) < thresholds(pthr)) ./ iterations;
    
    fix_fpr(pthr) =  sum(SUM.fixpvals_fpr(:,2) < thresholds(pthr)) ./ iterations;
    fix_pow(pthr) =  sum(SUM.fixpvals_pow(:,2) < thresholds(pthr)) ./ iterations;
    
    rand_fpr(pthr) =  sum(SUM.randpvals_fpr(:,2) < thresholds(pthr)) ./ iterations;
    rand_pow(pthr) =  sum(SUM.randpvals_pow(:,2) < thresholds(pthr)) ./ iterations;
    
    LRTrand_fpr(pthr) =  sum(SUM.LRTrandpvals_fpr(:,2) < thresholds(pthr)) ./ iterations;
    LRTrand_pow(pthr) =  sum(SUM.LRTrandpvals_pow(:,2) < thresholds(pthr)) ./ iterations;
       
    covfix_fpr(pthr) =  sum(SUM.fixpvals_fpr(:,3) < thresholds(pthr)) ./ iterations;
    covfix_pow(pthr) =  sum(SUM.fixpvals_pow(:,2) < thresholds(pthr)) ./ iterations;
    
end

hold on;
clear h

scatter(fix_fpr, fix_pow, 34, log(thresholds), 'filled');
scatter(covfix_fpr, covfix_pow, 34, log(thresholds), 'filled');
scatter(LRTrand_fpr, LRTrand_pow, 34, log(thresholds), 'filled');

h(1) = plot(sort(fix_fpr), sort(fix_pow), 'k');
h(2) = plot(sort(covfix_fpr), sort(covfix_pow), 'k-.');
h(3) = plot(sort(LRTrand_fpr), sort(LRTrand_pow), 'k:');


wh = find(thresholds - .05 == min(abs(thresholds - .05)));
plot([fix_fpr(wh) fix_fpr(wh)], [0 fix_pow(wh)], 'k--');
text(fix_fpr(wh) + .02, .1, 'p < .05','FontSize', 18);

wh = find(thresholds - .001 == min(abs(thresholds - .001)));
plot([fix_fpr(wh) fix_fpr(wh)], [0 fix_pow(wh)], 'k--');
text(fix_fpr(wh) + .01, .1, 'p < .001','FontSize', 18);

ylabel('Power'); xlabel('Est. Actual false positive rate')
legend(h, {'1st level fixed effects' '2nd level fixed effects' '2nd level random effects'});
set(gca, 'XLim', [-.01 .1]);

colorbar('vert');

