function [theta,HH,C,P,hrf,fit] = Anneal_Logit_allstim(theta0,t,tc,Run)
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
% further edited by Christian Waugh to include multiple trialtypes

numstim = length(Run);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initial values

iter = 3000*numstim;                               % Number of iterations
theta = theta0;                             % Set initial value for the parameter vector
h0 = cost_allstim(theta0,t,tc,Run);                 % Calculate cost of initial estimate
%LB = [0.05, 0, 0, 0.05, 2.5, 0.05, 5];      % Previous Lower bounds for parameters
%UB = [5, 5, 2, 2, 7.5, 2, 10];     
LB = [0.05, 1, 0, 0.01, 2.5, 0.05, 3];      % Lower bounds for parameters
UB = [6, 5, 2, 2, 7.5, 3, 10];           % Upper bounds for parameters
LB = repmat(LB, 1, numstim);
UB = repmat(UB, 1, numstim);

%
% These values may need tweaking depending on the individual situation.
% 

r1= 0.1;                                 % delta parameters
r1b= 0.1;                                % delta parameters
r2 = .5;                                 % T parameters
r3 = 0.1;                                % alpha parameters

t1 = [1 4];
t1b = [6];
t2 = [2 5 7];
t3 = [3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = zeros(1,7);
u = repmat(u, 1, numstim);

HH = zeros(1+iter,7*numstim);               % Keep track of theta_i
HH(1,:) = theta0;
P = zeros(1+iter,1);
C = zeros(1+iter,1);                % Keep track of the cost function
C(1) = h0;

cnt = 0;
for i=1:iter,
    
    T = 200/log(1+i);    %Temperature function (may require tweaking)
    th = zeros(1,7);
    th = repmat(th, 1, numstim);
    ind = 0;
    %time = zeros(numstim, 1);

    % Choose a new candidate solution theta_{i+1}, based on a random perturbation of the current solution of theta_{i}.
    % Check new parameters are within accepted bounds
    while ( (sum((LB-th)>0) + sum((th-UB)>0)) > 0) %|| (min(time) == 0),

        % Perturb solution
	
	for g = 0:(numstim-1)
        u(t1+g*7) = normrnd(0,r1,1,2);
        u(t1b+g*7) = normrnd(0,r1b,1,1);
        u(t2+g*7) = normrnd(0,r2,1,3);     
        u(t3+g*7) = normrnd(0,r3,1,1);
    end 
               
             
        % Update solution
        th = theta + u;
        ind = ind + 1;
        
        %include below if you want the times of inflection of the IL curves
        %to be in a specific order (i.e. T1 < T2 < T3)
        %for g = 1:numstim
          %  if th(g*7-5) < th(g*7-2) && th(g*7-2) < th(g*7)
           %     time(g) = 1;
           % else
            %    time(g) = 0;
            %end
        %end
        
        if(ind > 500), 
            warning('stuck!'); 
            th = theta; 
            %return; 
        end; 
    end;

    h = cost_allstim(th,t,tc,Run);
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
    hrf = zeros(length(t),numstim);
for g = 1:numstim
    hrf(:,g) = Get_Logit(theta(g*7-6:g*7),t);                   % Calculate HRF estimate (fit, given theta)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolve HRF and stick function
if nargout > 5
    len = length(Run{1});
    fitt = zeros(len,numstim);
for g = 1:numstim
    fits(:,g) = conv(Run{g}, hrf(:,g));
    fitt(:,g) = fits(1:len,g);
end
    fit = sum(fitt,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return
