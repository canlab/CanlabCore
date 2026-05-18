clear

% #########################################################################
% To validate the Eq. 2 of the paper and ensure it is identical the
% vectorised version is indetical to the matrix form (i.e. Eq. 1)
% SA, 2021
% #########################################################################

T = 500; 
nn = 2;
n1 = 1; n2 = 2;

% #########################################################################
% Simulate two time series under different conditions ---------------------
% #########################################################################

% Uncomment one of the following sections...

% uncorrelated white noise ################################################
%ts = randn(T,nn)';

% correlated white ########################################################
% If you choose to simulate this scenario, you need xDF added to your path
%
%rhosim = 0.7; 
%ts = corrautocorr([0 0],rhosim,eye(T),T);

% correlated AR1 ##########################################################
% If you choose to simulate this scenario, you need xDF added to your path
%
rhosim = 0.7; 
ARsim  = 0.5; 
ts = corrautocorr([0 0],rhosim,MakeMeCovMat(ARsim,T),T);

ts = ts-mean(ts,2); 


% xDF.m  ------------------------------------------------------------------
AC_ts=AC_fft(ts,T);

ac   = AC_ts (:,2:end-1);
ac_x = AC_ts(1,1:T-1);
ac_y = AC_ts(2,1:T-1);


xcf = xC_fft(ts,T);
xc_n      = flip(xcf(:,:,2:T-1),3);
xc_p      = xcf(:,:,T+1:end-1);

rho = corr(ts');

nLg    = T-2;        %if lag0 was the 0th element. Also, the ACF has T-1 dof.  
wgt     = (nLg:-1:1);
wgtm3   = reshape(repmat((repmat(wgt,[nn,1])),[nn,1]),[nn,nn,numel(wgt)]); 
wgtm2   = repmat(wgt,[nn,1]);

% #########################################################################
% Eq 2 --------------------------------------------------------------------
% #########################################################################

 FAST_New = ((T-1)*(1-rho.^2).^2 ...
     +   rho.^2 .* sum(wgtm3 .* (SumMat(ac.^2,nLg)  +  xc_p.^2 + xc_n.^2),3)...         %1 2 4
     -   2.*rho .* sum(wgtm3 .* (SumMat(ac,nLg)    .* (xc_p    + xc_n))  ,3)...         % 5 6 7 8
     +   2      .* sum(wgtm3 .* (ProdMat(ac,nLg)    + (xc_p   .* xc_n))  ,3))./(T^2);   % 3 9 
 %FAST_New(n1,n2)
   

[xcf_m,lags]  = crosscorr(ts(1,:),ts(2,:),T-1); %demean the time series
acx_n = fliplr(xcf_m(2:T));
acx_p = xcf_m(T:end-1);

%AACC=AC_fft(ts,T);

Sigma_X  = toeplitz(ac_x);
Sigma_Y  = toeplitz(ac_y);
Sigma_XY = (triu(toeplitz(acx_n))+tril(toeplitz(acx_p),-1));
Sigma_YX = Sigma_XY';

rho_m = rho(1,2);

% #########################################################################
% Eq 1 --------------------------------------------------------------------
% #########################################################################

%  ASAt = ((rho_m.^2./2) .* (trace(Sigma_X^2) + trace(Sigma_Y^2))... %1 & 2
%         + trace(Sigma_X*Sigma_Y) + rho_m^2 .* trace(Sigma_XY*Sigma_YX)... %3&4
%         - rho_m .* trace(Sigma_X*Sigma_XY) - rho_m .* trace(Sigma_X*Sigma_YX)... % 5&6
%         - rho_m .* trace(Sigma_Y*Sigma_XY) - rho_m .* trace(Sigma_Y*Sigma_YX)...%7&8
%         + trace(Sigma_XY^2))./T^2; %9
%num2str(ASAt) 

% Exact order as is in Eq 1 paper. 
ASAt0 = (T^(-2)).*( (rho_m.^2./2) .* trace(Sigma_X*Sigma_X) + (rho_m.^2./2) .* trace(Sigma_Y*Sigma_Y) ...         % 1 & 2
    + rho_m.^2 .* trace(Sigma_YX*Sigma_XY) + trace(Sigma_XY*Sigma_XY) + trace(Sigma_X*Sigma_Y) ... % 3&4
    - rho_m .* trace(Sigma_X*Sigma_XY) - rho_m .* trace(Sigma_YX*Sigma_X) - rho_m .* trace(Sigma_XY*Sigma_Y) - rho_m .* trace(Sigma_Y*Sigma_YX) ); %9   

disp(['Eq(1): ' num2str(ASAt0) ', Eq(2): ' num2str(FAST_New(1,2)) ])


%--------------------------------------------------------------------------
% Util function:
%--------------------------------------------------------------------------
function [SM0] = SumMat(Y0,T)
    %hopefully faster than Matlab's dumb sum of sum!
    %SA, Ox, 2018

    %becareful with this function, this is really tricky to use!
    if ~sum(ismember(size(Y0),T)); error('There is something wrong, mate!'); end
    if size(Y0,1) ~= T; Y0 = Y0'; end

    nn  = size(Y0,2);
    idx = find(triu(ones(nn),1))';
    SM0 = zeros(nn,nn,T);
    for i=idx
        [x,y]      = ind2sub([nn 1],i);
        SM0(x,y,:) = (Y0(:,x)+Y0(:,y));
        SM0(y,x,:) = (Y0(:,y)+Y0(:,x));
    end
end

function [SM0] = ProdMat(Y0,T)
    %hopefully faster than Matlab's dumb prod of prod!
    %SA, Ox, 2018

    %becareful with this function, this is really tricky to use!
    if ~sum(ismember(size(Y0),T)); error('There is something wrong, mate!'); end
    if size(Y0,1) ~= T; Y0 = Y0'; end

    nn  = size(Y0,2);
    idx = find(triu(ones(nn),1))';
    SM0 = zeros(nn,nn,T);
    for i=idx
        [x,y]      = ind2sub([nn 1],i);
        SM0(x,y,:) = (Y0(:,x).*Y0(:,y));
        SM0(y,x,:) = (Y0(:,y).*Y0(:,x));
    end
end

