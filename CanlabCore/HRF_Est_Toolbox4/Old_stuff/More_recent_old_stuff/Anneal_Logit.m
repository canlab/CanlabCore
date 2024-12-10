function [theta,HH,C,P,hrf,fit,e,param] = Anneal_Logit(theta0,t,tc,Run)
%
% [theta,HH,C,P] = Anneal_Logit(theta0,t,tc,Run)
%
% Estimate inverse logit (IL) HRF model using Simulated Annealing
% Creates fitted curve - 3 logistic functions to be summed together - from parameter estimates
%
% INPUT: theta0, t, tc, Run
% Run = stick function
% tc = time course
% t = vector of time points
% theta0 = initial value for the parameter vector
%
% By Martin Lindquist and Tor Wager
% Edited 12/12/06
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initial values

iter = 15000;                               % Number of iterations
theta = theta0;                             % Set initial value for the parameter vector
h0 = cost(theta0,t,tc,Run);                 % Calculate cost of initial estimate
LB = [0.05, 1, 0, 0.05, 5, 0, 10];      % Lower bounds for parameters
UB = [10, 15, 5, 10, 15, 5, 30];           % Upper bounds for parameters

%
% These values may need tweaking depending on the individual situation.
% 

r1= 0.001;                                 % A parameters
r1b= 0.001;                                % A parameters
r2 = 0.05;                                 % T parameters
r3 = 0.001;                                % delta parameters

t1 = [1 4];
t1b = [6];
t2 = [2 5 7];
t3 = [3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = zeros(1,7);

HH = zeros(1+iter,7);               % Keep track of theta_i
HH(1,:) = theta0;
P = zeros(1+iter,1);
C = zeros(1+iter,1);                % Keep track of the cost function
C(1) = h0;

cnt = 0;
for i=1:iter,
    
    T = 100/log(1+i);    %Temperature function (may require tweaking)
    th = zeros(1,7);
    ind = 0;

    % Choose a new candidate solution theta_{i+1}, based on a random perturbation of the current solution of theta_{i}.
    % Check new parameters are within accepted bounds
    while ( (sum((LB-th)>0) + sum((th-UB)>0)) > 0),

        % Perturb solution

        u(t1) = normrnd(0,r1,1,2);
        u(t1b) = normrnd(0,r1b,1,1);
        u(t2) = normrnd(0,r2,1,3);     
        u(t3) = normrnd(0,r3,1,1);

        % Update solution
        th = theta + u;
        ind = ind + 1;
        
        if(ind > 500), 
            warning('stuck!'); 
            return; 
        end; 
    end;

    h = cost(th,t,tc,Run);
    C(i+1) = h;
    delta = h - h0;
    
    % Determine whether to update the parameter vector.
    if (unifrnd(0,1) < min(exp(-delta/T),1)), 
        theta = th;
        h0=h;    
        cnt = cnt+1;
    end;

    HH(i+1,:) = theta;
    P(i+1) = min(exp(-delta/T),1);

end;

%cnt/iter

[a,b] = min(C);
theta = HH(b,:);
%h


% Additional outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get HRF for final model
if nargout > 4
    hrf = Get_Logit(theta(1:7),t);                   % Calculate HRF estimate (fit, given theta)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolve HRF and stick function
if nargout > 5
    len = length(Run);
    fit = conv(Run, hrf);
    fit = fit(1:len);
    e = tc - fit;
end

if nargout > 7
    [param] = get_parameters_logit(hrf,t,theta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return
