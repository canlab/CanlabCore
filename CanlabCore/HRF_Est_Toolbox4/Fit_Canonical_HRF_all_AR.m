function [HRF, Fit, Resid, Parameters, Results, df] = Fit_Canonical_HRF_all_AR(TC, TR, Run, T, p, Z, ar_lag)
% function [hrf, fit, e, param, info] = Fit_Canonical_HRF(tc,TR,Runs,T,p)
%
% Fits GLM using canonical hrf (with option of using time and dispersion derivatives)';
%
% INPUTS:
% 
% TC    - time course (voxel x time course)
% TR    - time resolution
% Runs  - expermental design
% T     - length of estimated HRF
% p     - Model type
% Z     - Nuisance covariates (time course x number of components)
% ar_lag - Number of lags in the AR model for residuals
%
% Options: p=1 - only canonical HRF
%          p=2 - canonical + temporal derivative
%          p=3 - canonical + time and dispersion derivative
% 
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual time course
% param - estimated amplitude, height and width
% info  - struct containing design matrices, beta values etc
%
% Created by Martin Lindquist on 07/11/24
% Last edited: 11/29/24

if nargin < 7
    ar_lag = 0; % Default value
end


numstim = length(Run);
len = length(Run{1});
t=1:TR:T;
num_vox = size(TC,1);
num_time = size(TC,2);

Resid = zeros(num_vox, num_time);
Fit = zeros(num_vox, num_time);
Results = cell(num_vox,1);
Parameters = zeros(num_vox, 3, numstim);
HRF = zeros(num_vox, length(t)+1, numstim);

X = zeros(len,p*numstim);
param = zeros(3,numstim);

[h, dh, dh2] = CanonicalBasisSet(TR);

for i=1:numstim
    v = conv(Run{i},h);
    X(:,(i-1)*p+1) = v(1:len);

    if (p>1)
        v = conv(Run{i},dh);
        X(:,(i-1)*p+2) = v(1:len);
    end

    if (p>2)
        v = conv(Run{i},dh2);
        X(:,(i-1)*p+3) = v(1:len);
    end
end
    
d = size(X,2);
X = [X Z (zeros(len,1)+1)];

df = size(X,2);

for v=1:num_vox

    y = TC(v,:)';
    [b, ~, ~, ~, ~, ~, ~, ~, df, ~, ~, ~] = fit_gls2(y, X, [], ar_lag);
    e = y - X * b; 
    fit = X(:, 1:d) * b(1:d);
    
    b = reshape(b(1:(p*numstim)),p,numstim)';
    bc = zeros(numstim,1);
    
    for i=1:numstim
        if (p == 1)
            bc(i) = b(i,1);
            H = h;
        elseif (p==2)
            bc(i) = sign(b(i,1))*sqrt((b(i,1))^2 + (b(i,2))^2); 
            H = [h dh];
        elseif (p>2)
            bc(i) = sign(b(i,1))*sqrt((b(i,1))^2 + (b(i,2))^2 + (b(i,3))^2);
            H = [h dh dh2];
        end    
    
    end
    
    hrf = H*b';
    
    for i=1:numstim
        param(:,i) = get_parameters2(hrf(:,i),1:length(t));
    end
    
    info ={};
    info.b = b;
    info.bc = bc;
    info.X = X;
    info.H =H;


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunctions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h, dh, dh2] = CanonicalBasisSet(TR)

len = round(30/TR);
xBF.dt = TR;
xBF.length= len;
xBF.name = 'hrf (with time and dispersion derivatives)';
xBF = spm_get_bf(xBF);

v1 = xBF.bf(1:len,1);
v2 = xBF.bf(1:len,2);
v3 = xBF.bf(1:len,3);

h = v1;
dh =  v2 - (v2'*v1/norm(v1)^2).*v1;
dh2 =  v3 - (v3'*v1/norm(v1)^2).*v1 - (v3'*dh/norm(dh)^2).*dh;

h = h./max(h);
dh = dh./max(dh);
dh2 = dh2./max(dh2);

end