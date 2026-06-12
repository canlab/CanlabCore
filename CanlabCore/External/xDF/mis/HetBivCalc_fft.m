function [CF,V,EDF,xAC,acpvals]=HetBivCalc_fft(Y,L,varargin)
%[CF,EDF,xAC]=HetBivCalc_fft(Y,L,varargin)
%
%   Super fast full-lag Bartlett's Correction Factor (BCF) calculation of 
%   multi-dimention matrices. 
%
%%%INPUTS:
%   Y: IxT or TxI BOLD voxel-wise/pacellated time series 
%   L: Length of the time series! Just to check the dimensions are all
%   sorted!
%
%   Optional Inputs:
%       Method: 'BHCF', 'CF0', 'CFx'. [default: BHCF]
%       lag: curbing on cross-covariance structure. [default: 1%]
%
%
%%%OUTPUTS:
%   BCF:  IxI matrix. Correction Factor selected as Method. \hat{N}=N/BCF
%   EDF:  IxI matrix. Effective Degree of Freedom.
%   xAC:  Full-lag Autocorrelation. IxT-1, Becareful about memory! 
%   NB! if the correction factor is below 1, then it is forced to 1.
%%%DEPENDENCIES: 
%   AC_fft.m 
%   xC_fft.m
%
%%%REFERENCES:
%
%   Afyouni & Nichols, 2017
%   Bayley & Hammersley, 1946
%   Richardson & Clifford, 1991
%_________________________________________________________________________
% Soroosh Afyouni, NISOx.org, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

restricttheones=1; trimflag = 0; 
%change this to 1 if you want to restrict them to min 1.
%78 85 98 these lines make sure you never exceed the original df. They are
%commented atm, return them back after the simulation tests. 

%This should be approximately equivalent to CF=HetBiv_fft(t,N)
%trace(toeplitz(AC_fft(t(1,:),N))*toeplitz(AC_fft(t(2,:),N)))./N

if size(Y,2)~=L
    Y=Y'; %IxT
end
I = size(Y,1);

if sum(strcmpi(varargin,'Method'))
   CFm = varargin{find(strcmpi(varargin,'Method'))+1};
    if strcmpi(CFm,'BHCF')
        CFmethod=1;
        lagcc=0;
    elseif strcmpi(CFm,'CF0')
        CFmethod=2;
        lagcc=0;
    elseif strcmpi(CFm,'CFx')
        CFmethod=3;
        if sum(strcmpi(varargin,'lag'))
           lagcc = varargin{find(strcmpi(varargin,'lag'))+1};
        else
           lagcc = round(L*0.01); %because we found that if you keep it short it would work better.
           if ~lagcc; lagcc=1; end;
        end
    else
        error('Unknown CF methods. Choose between BHCF, BC0, BCx');
    end
else
   CFmethod=1;
   lagcc=0;
end

if sum(strcmpi(varargin,'Trim'))
    Trim_str = varargin{find(strcmpi(varargin,'Trim'))+1};
    if strcmpi(Trim_str,'on')
        trimflag=1;
    elseif strcmpi(Trim_str,'off')
        trimflag=0;
    elseif strcmpi(Trim_str,'taper') 
        trimflag=2;
    else
        error('Trim should be followed by Tucky, taper or OFF.')
    end
end

%if sum(strcmpi(varargin,'Whiten'))
%    %disp('--WHITENED.')
%    Y = WhiteningCorr(Y,L,'Method','SQRTM');
%end

xAC      = AC_fft(Y,L);

acpvals=zeros(size(xAC));
if trimflag==1
%    disp('--TRIM ON.')
%     error('JUST DONT FOR NOW!!')
%     %----Detecting the sig AC lags: 
%     varacf   = (1+2.*sum(xAC(:,1:(L/5)).^2,2))./L; %From Anderson's p8: variance of a.c.f        
%     zs       = xAC./sqrt(varacf);     %z-scores
%     acpvals0 = 2.*normcdf(-abs(zs));  %pvals
%     %FDR
%     for i=1:I; acpvals(i,:) = fdr_bh(acpvals0(i,:)); end; %FDR
%     xAC     = acpvals.*xAC;%filter the AC function

M              = round(varargin{find(strcmpi(varargin,'on'))+1});
acpvals(:,1:M) = 1;
xAC            = acpvals.*xAC;

elseif trimflag == 0
    %disp('--TRIM OFF.')
    acpvals = ones(size(xAC));
    xAC     = acpvals.*xAC;
elseif trimflag==2
    %disp('--TURKEY WINDOW TAPERED.')
    M              = round(varargin{find(strcmpi(varargin,'taper'))+1});
    %M             = round(sqrt(L));
    xAC_tmp        = zeros(size(xAC));
    xAC_tmp(:,1)   = 1;  % becareful here! you'll remove the lag-0 further below!
    xAC_tmp(:,2:M) = (1+cos([2:M].*pi./M))./2.*xAC(:,2:M);
    
    xAC            = xAC_tmp; 
    clear *_tmp acpvals;
end

xAC(:,1) = []; %because we later take care of that little lag-0!

%----
nLg    = L-1;
wgt    = (nLg:-1:1);
CF     = wgt.*xAC*xAC'; %pfff
CovSxy = (L+2*(CF));
CF     = CovSxy./L;
V      = CovSxy./L.^2;
% this dude above should be trace(toeplitz(ac100)*toeplitz(ac200))./N

if CFmethod==1
    if restricttheones; CF(CF<1) = 1; end; %ensure that there are no node with BCF smaller than 1.
    EDF = L./CF; 
    return; 
end

if CFmethod==2
    rho    = corr(Y');    
    CovSxy = CovSxy+(L*rho.^2);
    CF     = (CovSxy)./L;
    V      = (CovSxy)./L.^2;    
    if restricttheones; CF(CF<1) = 1; end; %ensure that there are no node with BCF smaller than 1.
    EDF    = L./CF; 
end

if CFmethod==3
    xC  = xC_fft(Y,L,'lag',lagcc); %IxIxlagcc
    wgt = (nLg-abs(-(lagcc-1):lagcc-1)); %1xlagcc
    %size(wgt), size(xC)
    xC  = bsxfun(@times,xC.^2,reshape(wgt,1,1,numel(wgt))); %here, wgt is multiplied by each plain of IxI of xC!
    xC  = sum(xC,3);
    xC  = xC+triu(xC,1)';
    CF  = CF+xC./L;
    
    if restricttheones; CF(CF<1) = 1; end; %ensure that there are no node with BCF smaller than 1.
    
    EDF = L./CF;
end 


