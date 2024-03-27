function [hrf, fit, e, param, info] = Fit_Canonical_HRF(tc, TR, Run, T, p, varargin)
% function [hrf, fit, e, param, info] = Fit_Canonical_HRF(tc,TR,Runs,T,p)
%
% Fits GLM using canonical hrf (with option of using time and dispersion derivatives)';
%
% INPUTS:
% 
% tc    - time course
% TR    - time resolution
% Runs  - expermental design
% T     - length of estimated HRF
% p     - Model type
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
% Created by Martin Lindquist on 10/02/09
% Last edited: 05/17/13 (ML)

% Edited to allow passing in custom design matrix e.g., from SPM. Will
% assume that the first regressor columns of the design matrix pertain to
% the regressors in Run: Michael Sun, Ph.D. 02/20/2024

[h, dh, dh2] = CanonicalBasisSet(TR);
%tc = tc';
d = length(Run);
len = length(Run{1});
% Generate a design matrix
t=1:TR:T;

% Import your own design matrix
if ~isempty(varargin)
    X=varargin{1};
    if numel(varargin)>1
        for i = 1:numel(varargin)
            if strcmpi(varargin{i}, 'invertedDX')
                PX=varargin{i+1};
                b=PX*tc;
            end 
        end
    end
    if ~exist('b', 'var')
        b = pinv(X)*tc;
    end

    e = tc-X*b;
    fit = X*b;
    
    % Be careful here. if p>1, make sure Run includes derivatives so there
    % are p*task regressors.
    b = reshape(b(1:numel(Run)),p,d)'; % Extract my own regressors
    bc = zeros(d,1);

else
   
    % Constructing the Design Matrix X:
    X = zeros(len,p*d);
    param = zeros(3,d);
    
    for i=1:d,
        v = conv(Run{i},h);
        X(:,(i-1)*p+1) = v(1:len);
    
        % Computing the first derivative
        if (p>1)
            v = conv(Run{i},dh);
            X(:,(i-1)*p+2) = v(1:len);
        end
    
        % Computing the second derivative
        if (p>2)
            v = conv(Run{i},dh2);
            X(:,(i-1)*p+3) = v(1:len);
        end
    end
    
    % This line adds an intercept
    X = [(zeros(len,1)+1) X];
    PX = pinv(X);
    b = PX*tc;
    e = tc-X*b;
    fit = X*b;
   
    b = reshape(b(2:end),p,d)';
    bc = zeros(d,1);
end

for i=1:d,
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

for i=1:d,
    param(:,i) = get_parameters2(hrf(:,i),1:length(t));
end;


info ={};
info.b = b;
info.bc = bc;
info.DX = X;
info.PX = PX;
info.H =H;

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