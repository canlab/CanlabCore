%function [hrf, fit, e, param] = Fit_sFIR(tc, TR, Run, T, mode)
function [HRF, Fit, Resid, Parameters, Results, df] = Fit_sFIR_all(TC, TR, Runc, T, mode, Z)
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
% Created by Martin Lindquist on 07/11/24
% Last edited: 

numstim = length(Runc);
len = length(Runc{1});
t=1:TR:T;
tlen = length(t);
num_vox = size(TC,1);
num_time = size(TC,2);

Resid = zeros(num_vox, num_time);
Fit = zeros(num_vox, num_time);
Results = cell(num_vox,1);
Parameters = zeros(num_vox, 3, numstim);
HRF = zeros(num_vox, tlen, numstim);

Runs = zeros(len,numstim);
for i=1:numstim
    Runs(:,i) = Runc{i};
end

[DX] = tor_make_deconv_mtx3(Runs,tlen,1);
X = [DX Z];

df = size(X,2);


if mode == 1

    C=(1:tlen)'*(ones(1,tlen));
    h = sqrt(1/(7/TR));                       % 7 seconds smoothing - ref. Goutte

    v = 0.1;
    sig = 1;

    R = v*exp(-h/2*(C-C').^2);
    RI = inv(R);
    MRI = zeros(numstim*tlen+1);
    for i=1:numstim
        MRI(((i-1)*tlen+1):(i*tlen),((i-1)*tlen+1):(i*tlen)) = RI;
    end


    f = inv(blkdiag(DX'*DX+sig^2*MRI, Z'*Z))*X';
   
    for v=1:num_vox
        y = TC(v,:)';
        b = f*y;
        fit = X*b;
        e = y - fit; 

        hrf =zeros(tlen,numstim);
        param = zeros(3,numstim);

        for i=1:numstim
            hrf(:,i) = b(((i-1)*tlen+1):(i*tlen))';
            param(:,i) = get_parameters2(hrf(:,i),(1:tlen));
        end

        info ={};
        Resid(v,:) = e';
        Fit(v,:) = fit';
        Results{v} = info;
        Parameters(v,:,:) = param;
        HRF(v,:,:) = hrf;

        disp(v);

    end

elseif mode == 0

    pX = pinv(X);
    for v=1:num_vox
        y = TC(v,:)';
        b = pX*y;
        fit = X*b;
        e = y - fit; 
        
        hrf =zeros(tlen,numstim);
        param = zeros(3,numstim);

        for i=1:numstim
            hrf(:,i) = b(((i-1)*tlen+1):(i*tlen))';
            param(:,i) = get_parameters2(hrf(:,i),(1:tlen));
        end

        info ={};
        Resid(v,:) = e';
        Fit(v,:) = fit';
        Results{v} = info;
        Parameters(v,:,:) = param;
        HRF(v,:,:) = hrf;

        disp(v)
    end    
end



end