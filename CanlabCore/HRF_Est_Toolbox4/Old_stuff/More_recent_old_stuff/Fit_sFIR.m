function [hrf, fit, e, param] = Fit_sFIR(tc,TR,Run,T,mode)
% function [hrf, fit, e, param] = Fit_sFIR(tc,TR,Runs,T,mode)
%
% Fits FIR and smooth FIR model  
%
% INPUTS:
% 
% tc    - time course
% TR    - time resolution
% Runs  - expermental design
% T     - length of estimated HRF
% mode  - FIR or smooth FIR
%   options:
%       0 - standard FIR 
%       1 - smooth FIR
% 
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual time course
% param - estimated amplitude, height and width
%
% Created by Martin Lindquist on 10/02/09
% Last edited: 05/26/10 (ML)

numstim = length(Run);
len = length(Run{1});


Runs = zeros(len,numstim);
for i=1:numstim,
    Runs(:,i) = Run{i};
end;

[DX] = tor_make_deconv_mtx3(Runs,T,1);

% DX2 = DX(:,1:T); 
% num = T;

if mode == 1

    C=(1:T)'*(ones(1,T));
    h = sqrt(1/(7/TR));                       % 7 seconds smoothing - ref. Goutte

    v = 0.1;
    sig = 1;

    R = v*exp(-h/2*(C-C').^2);
    RI = inv(R);
    MRI = zeros(numstim*T+1);
    for i=1:numstim,
        MRI(((i-1)*T+1):(i*T),((i-1)*T+1):(i*T)) = RI;
    end;

    b = inv(DX'*DX+sig^2*MRI)*DX'*tc;
    fit = DX*b;
    e = tc - DX*b; 

elseif mode == 0

    b = pinv(DX)*tc;
    fit = DX*b;
    e = tc - DX*b;
    
end


hrf =zeros(T,numstim);
param = zeros(3,numstim);

for i=1:numstim,
    hrf(:,i) = b(((i-1)*T+1):(i*T))';
    param(:,i) = get_parameters2(hrf(:,i),T);
end;

end