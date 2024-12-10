%function [hrf, fit, e, param] = Fit_Spline_all(tc, TR, Run, T)
function [HRF, Fit, Resid, Parameters, Results, df] = Fit_Spline_all(TC, TR, Run, T, Z)

% function [hrf, fit, e, param] = Fit_Spline(tc, TR, Run, T)
%
% Fits Spline model  
%
% INPUTS:
% 
% TC    - time course (voxel x time course)
% TR    - time resolution
% Runs  - expermental design
% T     - length of estimated HRF
% Z     - Nuisance covariates (time course x number of components)
% 
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual time course
% param - estimated amplitude, height and width
%
% Created by Martin Lindquist on 07/11/24
% Last edited: 

numstim = length(Run);  % Number conditions
len = length(Run{1});   % length of run
t=1:TR:T;                   
tlen = length(t);       % Number of time points in HRF
num_vox = size(TC,1);
num_time = size(TC,2);

Resid = zeros(num_vox, num_time);
Fit = zeros(num_vox, num_time);
Results = cell(num_vox,1);
Parameters = zeros(num_vox, 3, numstim);
HRF = zeros(num_vox, tlen, numstim);

param = zeros(3,numstim);



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

X = [Wi Z (zeros(len,1)+1)];      

df = size(X,2);

pinvX = pinv(X);

% Fit model
for v=1:num_vox

    y = TC(v,:)';
    b = pinvX*y;
    e = y-X*b;
    fit = X*b;
    
    b = reshape(b(1:(K*numstim)),K,numstim);

    % Get parameters
    
    hrf = full(B)*b;
    
    for i=1:numstim
        param(:,i) = get_parameters2(hrf(:,i),1:tlen);
    end

    info ={};
    info.b = b;
    info.X = X;
    info.H = B;

    Resid(v,:) = e';
    Fit(v,:) = fit';
    Results{v} = info;
    Parameters(v,:,:) = param;
    HRF(v,:,:) = hrf;

    disp(v)

end







end

% END MAIN FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
