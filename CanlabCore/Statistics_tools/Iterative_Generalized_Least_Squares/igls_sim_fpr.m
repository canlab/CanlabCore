%s = rand('twister');
%rand('twister', s);

len = 200;              % number of 1st level observations (time points within subjects)
sub = 20;               % number of second-level units (subjects)

fixed_slope = 0.5;      % population fixed-effect slope
fixed_int = 0.3;        % pop. fixed-effect intercept

rand_slope = 0.5;         % 'true' random effect standard dev: slope
rand_int = 0.5;           % 'true' random effect standard dev: intercept

noise_sigma = 0.5;      % measurement error/unexplained std of 1st level observations
between_sigma = 0.5;    % measurement error/unexplained std of 2nd level covariate

% Generate simulated data
% 1) create predictors
x = zeros(len,sub);
y = x;
x(11:20,:) = 1;                   
x(111:120,:) = 1;                  

% 2) Create instances of effects for each subject
c = normrnd(fixed_slope, rand_slope, sub, 1);       % slope between-subjects variations
d = normrnd(fixed_int, rand_int, sub, 1);                  % intercept between-subjects variations

% 3) Create a between-ss covariate that is related to the slope but not the
% intercept
covt = scale(c + normrnd(0, between_sigma, sub, 1), 1);

% 4) Create a between-ss covariate that is related to the intercept but not the
% intercept
covti = scale(d + normrnd(0, between_sigma, sub, 1), 1);

corrcoef([c d covt covti])

% Add between-subjects error (random effects) and measurement noise
% (within-subjects error)

for i=1:sub
    y(:,i) = d(i) + c(i) .* x(:,i) + normrnd(0, noise_sigma, len, 1);
end

%% Run the model - without a 2nd-level covariate
out = igls(y, x);  % for igls
fprintf('\t\t%s\t%s\t%s\t\n', 'Effect', 'Intercept', 'Slope');
disp('True random-effect variances:'); 
fprintf('\t\t\t%3.3f\t%3.3f\n', [rand_slope rand_int]);

disp('Input random-effect variances: '); 
fprintf('\t\t\t%3.3f\t%3.3f\n', std([d c]))

disp('Est.  random-effect variances: '); 
fprintf('\t\t\t%3.3f\t%3.3f\n', sqrt(out.betastar)');

%% Run the model - with an unrelated 2nd-level covariate
out = igls(y, x, 'covariate', (1:20)');  % for igls
fprintf('\t\t%s\t%s\t%s\t\n', 'Effect', 'Intercept', 'Slope');
disp('True random-effect variances:'); 
fprintf('\t\t\t%3.3f\t%3.3f\n', [rand_slope rand_int]);

disp('Input random-effect variances: '); 
fprintf('\t\t\t%3.3f\t%3.3f\n', std([d c]))

disp('Est.  random-effect variances: '); 
fprintf('\t\t\t%3.3f\t%3.3f\n', sqrt(out.betastar)');

%% Run the model - with a 2nd-level covariate related to slope
out = igls(y, x, 'covariate', covt);  % for igls
fprintf('\t\t%s\t%s\t%s\t\n', 'Effect', 'Intercept', 'Slope');
disp('True random-effect variances:'); 
fprintf('\t\t\t%3.3f\t%3.3f\n', [rand_slope rand_int]);

disp('Input random-effect variances: '); 
fprintf('\t\t\t%3.3f\t%3.3f\n', std([d c]))

disp('Est.  random-effect variances: '); 
fprintf('\t\t\t%3.3f\t%3.3f\n', sqrt(out.betastar)');

%% Run the model - with a 2nd-level covariate related to intercept
out = igls(y, x, 'covariate', covti);  % for igls
fprintf('\t\t%s\t%s\t%s\t\n', 'Effect', 'Intercept', 'Slope');
disp('True random-effect variances:'); 
fprintf('\t\t\t%3.3f\t%3.3f\n', [rand_slope rand_int]);

disp('Input random-effect variances: '); 
fprintf('\t\t\t%3.3f\t%3.3f\n', std([d c]))

disp('Est.  random-effect variances: '); 
fprintf('\t\t\t%3.3f\t%3.3f\n', sqrt(out.betastar)');

%%
clear *pvals* *fpr*

mytype = 'i';
iterations = 200;

len = 100; sub = 10;

fixpvals = zeros(iterations, 2);
randpvals = zeros(iterations, 2);
LRTrandpvals = zeros(iterations, 2);
chic = zeros(iterations, 1);
chid = zeros(iterations, 1);

c_fixed = 0;  % slope pop. average
d_fixed = 0;   % intercept pop. average
c_rand = 0;    % slope std across Ss
d_rand = 0;    % intercept std.
noise_std = 2.0;  % within-subjects noise std

x = zeros(len,sub);
y = x;
x(1:2:end,:) = 1;                   % create signal

% --------------------------------------------------------------------------
% This block: No true random effects, so test false positive rate for
% p_randvariance

verbstr = 'noverbose';

for i = 1:iterations
  
    
    c = normrnd(c_fixed,c_rand,sub,1);       % slope between-subjects variations
    d = normrnd(d_fixed,d_rand,sub,1);         % intercept between-subjects variations
    
    % Add between-subjects error (random effects) and measurement noise
    % (within-subjects error)
           y = zeros(len, sub); 
    for s = 1:sub
        y(:,s) = d(s) + c(s) .* x(:,s) + normrnd(0,noise_std,len,1);
    end
    out = igls(y, x, 'type', mytype, verbstr);  % for igls
 
    fpr_converged(i) = out.isconverged;
    
    
    fpr_betastar(i, :) = out.betastar';
    
    fixpvals(i, :) = out.p;
    randpvals(i, :) = out.p_randvariance'; %[out.p_randvariance_d out.p_randvariance_c];
    
    LRTrandpvals(i, :) = out.pLRT_randvariance';
    
    if mod(i, 10) == 0, fprintf(1, '%3.0f ', i); end
    
    verbstr = 'noverbose';
end
fprintf('\n')

%pvals_fpr = pvals;
alph = .1;

fixfpr =  sum(fixpvals < alph) ./ iterations;
randfpr =  sum(randpvals < alph) ./ iterations;
LRTfpr =  sum(LRTrandpvals < alph) ./ iterations;

fprintf('\tfixed\t\trand fx\t\t\n');
fprintf('TPR/FPR @ alpha = %3.3f:\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f\n', alph, fixfpr, randfpr, LRTfpr);