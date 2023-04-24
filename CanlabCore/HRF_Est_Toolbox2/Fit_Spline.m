function [hrf, fit, e, param] = Fit_Spline(tc, TR, Run, T)
% function [hrf, fit, e, param] = Fit_Spline(tc, TR, Run, T)
%
% Fits Spline model  
%
% INPUTS:
% 
% tc    - time course
% TR    - time resolution
% Runs  - expermental design
% T     - length of estimated HRF
% 
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual time course
% param - estimated amplitude, height and width
%
% Created by Martin Lindquist on 03/06/23

numstim = length(Run);  % Number conditions
len = length(Run{1});   % length of run
t=1:TR:T;                   
tlen = length(t);       % Number of time points in HRF

K = 8;                         % Number of b-spline basis (This cxan be set to specification)
norder = 4;                    % Order of b-spline basis (This cxan be set to specification)

% Create design matrix
basis = create_bspline_basis([0,tlen], K+3, norder);    
B = eval_basis((1:tlen),basis);
B = B(:,3:end-1);

Wi = zeros(len, numstim*K);
for j=1:numstim  
    Wji = tor_make_deconv_mtx3(Run{j},tlen,1);    
    Wi(:,(j-1)*K+1:j*K) = Wji(:,1:tlen)*B;    
end

X = [ones(len,1) Wi];      

% Fit model
b = pinv(X)*tc;
e = tc-X*b;
fit = X*b;

b2 = reshape(b(2:end),K,numstim);


% Get parameters

hrf = B*b2;

for i=1:numstim
    param(:,i) = get_parameters2(hrf(:,i),1:length(t));
end

end

% END MAIN FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
