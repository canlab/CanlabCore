function [hrf, fit, e, param] = Fit_sFIR(tc,TR,Runs,T,mode)
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

[DX] = tor_make_deconv_mtx3(Runs,T,1);
DX2 = DX(:,1:T); 
num = T;

if mode == 1

    C=(1:num)'*(ones(1,num));
    h = sqrt(1/(7/TR));                       % 7 seconds smoothing - ref. Goutte

    v = 0.1;
    sig = 1;

    R = v*exp(-h/2*(C-C').^2);
    RI = inv(R);

    b = inv(DX2'*DX2+sig^2*RI)*DX2'*tc;
    fit = DX2*b;
    e = tc - DX2*b; 

elseif mode == 0

    b = pinv(DX)*tc;
    fit = DX*b;
    e = tc - DX*b;
    
end

hrf = b;
param = get_parameters2(b,T);

end