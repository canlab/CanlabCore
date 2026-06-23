function [cf,edf]=BCFCalc(x1,x2,varargin)
%[cf,dfn]=BCFCalc(x1,x2,Howfar)
%   
%   Calculates the effective degrees of freedom for correlation coefficient
%   of two give time series. 
%
%%%INPUTS
%   x1  : first time series as a vector
%   x2  : second time series as a vector
%%%OUTPUTS
%   cf  : is the correction factor
%   edf : is the effective degrees of freedom (i.e. N/cf)
%
%%%NOTES
%
%   1) time series, x1 and x2, should be mean-zeroed!
%   2) Pyper and Peterman recommend *not* to use Garrett-Petrie on short
%   time series. 
%   3) They also recommend Chelton method over Kope-Botsford. However they
%   seems to be quite similar in results.
%   4) Robust AC is NOT NOT NOT recommanded.
%
%%%REFERENCE
%
%   Pyper, B. J., & Peterman, R. M. (1998). Comparison of methods to 
%   account for autocorrelation in correlation analyses of fish data, 
%   2140, 2127?2140.
%
%   Chelton, D.B. 1984. Commentary: short-term climatic variability 
%   in the Northeast Pacific Ocean. In The influence of ocean conditions
%   on the production of salmonids in the North Pacific. 
%   Edited by W. Pearcy. Oregon State University Press, Corvallis, Oreg. pp.
%
%   Kope, R.G., and Botsford, L.W. 1990. Determination of factors affecting
%   recruitment of chinook salmon, Oncorhynchus tshawyt- scha, in central 
%   California. Fish. Bull. U.S. 88: 257?269.
%
%   Garrett, C., and Petrie, B. 1981. Dynamical aspects of the flow through
%   the Strait of Belle Isle. J. Phys. Oceanogr. 11: 376?393.
%
%
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________
%% Check Param

ndpr=length(x1);

if sum(strcmpi(varargin,'Howfar'))
   lagkey           =   varargin{find(strcmpi(varargin,'Howfar'))+1};
else
   lagkey   =    round(ndpr./5); % 20%  
end

if sum(strcmpi(varargin,'ACtype'))
    ACm           =   varargin{find(strcmpi(varargin,'ACtype'))+1};
    if strcmpi(ACm,'Unbiased')
        ACmfalg=1;
    elseif strcmpi(ACm,'robust')
        ACmfalg=2;
    else
        error('Unknown method!')
    end
else
    ACmfalg=0;
end

if sum(strcmpi(varargin,'Interaction'))
   InterMeth           =   varargin{find(strcmpi(varargin,'Interaction'))+1};
    if strcmpi(InterMeth,'GarretPetrie')
        InterMethflag=1;
    else
    error('Unknown method!')
    end
else
   InterMethflag   =    0; %continues with Pyper&Peterman
end

if sum(strcmpi(varargin,'Approximation'))
   Method           =   varargin{find(strcmpi(varargin,'Approximation'))+1};
    if strcmpi(Method,'Bartlett')
        wghtflag=0;
    else
        error('Unknown method!')
   end
else
    wghtflag=1; %continues with Chelton.
end

if size(x1,1)~=(ndpr)  
    disp(['TS1 timeseries transposed!'])
    x1=x1';
end
if size(x2,1)~=(ndpr)
    disp(['TS2 timeseries transposed!'])
    x2=x2';
end

%% 

if  ACmfalg==0 %Biased
    ac_x=autocorr(x1,lagkey);
    ac_y=autocorr(x2,lagkey);
    
elseif ACmfalg==1 %Unbiased
    ac_x=(ndpr./(ndpr-(0:lagkey)))'.*autocorr(x1,lagkey);
    ac_y=(ndpr./(ndpr-(0:lagkey)))'.*autocorr(x2,lagkey);
    
elseif ACmfalg==2 %Robust
    ac_x(1)=1; ac_y(1)=1;
    for jj=2:lagkey+1
        ac_x(jj,1)= madicc(x1(1:end-jj),x1(jj+1:end));
        ac_y(jj,1)= madicc(x2(1:end-jj),x2(jj+1:end));
    end
end

if InterMethflag
    rho_xy=cross(ac_x(2:end)',ac_y(2:end)');
elseif ~InterMethflag
    rho_xy=ac_x(2:end).*ac_y(2:end);
end

if wghtflag
    wght=(ndpr-(1:lagkey))./ndpr;
elseif ~wghtflag
    wght=1;
end

cf=1+2*sum(wght.*rho_xy');
edf=ndpr./cf;


function rmad = madicc(x,y)
% Median Absolute Deviation Intraclass Correlation Coefficient
%TEN, UoW, 2017

I=find(all(~isnan([x(:) y(:)]),2));
if isempty(I)
  rmad=NaN;
else
  mx    = median(x(I));
  my    = median(y(I));
  Sp    = (x(I)-mx) + (y(I)-my);
  Sm    = (x(I)-mx) - (y(I)-my);
  madSp = median(abs(Sp-median(Sp)));
  madSm = median(abs(Sm-median(Sm)));
  if madSp==0 && madSm==0
    rmad = NaN;
  else
    rmad = (madSp^2 - madSm^2)/(madSp^2 + madSm^2);
  end
end
    