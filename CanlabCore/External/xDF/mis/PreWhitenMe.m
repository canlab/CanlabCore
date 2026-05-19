function [wY,W]=PreWhitenMe(Y,T,varargin)
%   [wY,W_tmp]=PreWhitenMe(Y,T,varargin)
%   Prewhiten the time series; whiten the autocorrelation.   
%   
%   Y: can be a vector of size T or a matrix of size KxT comprised of K 
%   time series. Each time series is pre-whitened independently. 
%   T: number of datapoints; for sake of sanity check!
%   
%%%OUTPUTS:
%   wY : prewhitened time series 
%   W  : prewhitening matrices, if you are prewhitening multiple time
%   series it is gonna be a TxTxK 3D matrix (becareful with memory!)
%   
%%%EXAMPLE:
%  Use a single tukey tapering (Woolrich et al 2001) & Cholesky decomposition:
%  [wY,W_tmp]=PreWhitenMe(Y,T,'taper','tukey',sqrt(T),'DM','Cholesky');
%   
%  Use an AR2 model & SVD decomposition:
%  [wY,W_tmp]=PreWhitenMe(Y,T,'truncate',2,'DM','svd');
%
%%%NOTES:
%  Do not forget to specify your model via 'taper' or 'truncation' option.
%
%%%REFERENCES:
%  Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. 
% "Effective Degrees of Freedom of the Pearson's Correlation Coefficient 
%  under Serial Correlation." bioRxiv (2018): 453795.
%
%  Woolrich, M. W., Ripley, B. D., Brady, M., & Smith, S. M. (2001). 
%  Temporal autocorrelation in univariate linear modeling of FMRI data. 
%  NeuroImage, 14(6), 1370?86. http://doi.org/10.1006/nimg.2001.0931
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2018
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

%clear
% T = 1000;
% ts  = corrautocorr([0 0],0.9,0.4,T);
% Y   = ts(1,:);
% DM = 'SVD';

    if  size(Y,2) ~= T
        Y = Y';
        warning('Oi!')
    end
    
    Nts = size(Y,1);
    
    mY  = Y-mean(Y,2); 
    Yac = AC_fft(Y,T);
    
    nLg = T-2;
    
%-------READING
    if sum(strcmpi(varargin,'DM'))
       DM = lower(varargin{find(strcmpi(varargin,'DM'))+1});
    else
       DM = 'svd'; %because you don't need to be worry about PSD matrices. 
    end

    if sum(strcmpi(varargin,'taper'))
        mth = varargin{find(strcmpi(varargin,'taper'))+1};
        if strcmpi(mth,'tukey')
    %Tukey Tappering----------------------------------------
            M = round(varargin{find(strcmpi(varargin,'tukey'))+1});
            if isempty(M)
                M = sqrt(T); 
            end;
            %figure; plot(Yac(1,:)); hold on; plot(tukeytaperme(Yac(1,:),T,M));
            for in=1:Nts    
                Yac(in,:) = tukeytaperme(Yac(in,:),T,M);
            end
        else
            error('Umm, something is wrong in tapering options!')
        end
    %Shrinking------------------------------------------------
    elseif sum(strcmpi(varargin,'truncate'))
            mth = round(varargin{find(strcmpi(varargin,'truncate'))+1});
            if isnumeric(mth)
                for in=1:Nts    
                    Yac(in,:) = curbtaperme(Yac(in,:),T,mth);
                end
            elseif strcmpi(mth,'adaptive')
                for in=1:Nts    
                    Yac(in,:)= shrinkme(Yac(in,:),T);
                end
            else
                error('Umm, something is wrong in truncation options!')
            end
    %--------------------------------------------------------------------------
    else
        error('choose one of these; shrink | tukey | curb as tapering option.')
    end

%-----PREWHITENING

    for i=1:Nts
        iacm_tmp    = inv(toeplitz(Yac(i,:)));
        switch DM
            case 'cholesky'
                %disp(['Matrix decomp method: ' DM])
                W_tmp  = chol(iacm_tmp);
            case 'svd'
                %disp(['Matrix decomp method: ' DM])
                [U_tmp,S_tmp,~] = svd(iacm_tmp);
                W_tmp = U_tmp * sqrt(S_tmp) * U_tmp';
            case 'sqrtm'
                %disp(['Matrix decomp method: ' DM])
                error('Yeah mate, you dont have that much time!');
        end
        wY(i,:)  = W_tmp*mY(i,:)';
        W(:,:,i) = W_tmp; 
        clear *_tmp;
    end
end

%----OTHER FUNCTIONS

function srnkd_ts=shrinkme(acs,T)
%Shrinks the *early* bucnhes of autocorr coefficients beyond the CI.
%Yo! this should be transformed to the matrix form, those fors at the top
%are bleak!
%
%SA, Ox, 2018
    if ~sum(ismember(size(acs),T)); error('There is something wrong, mate!'); end
    if size(acs,2) ~= T; acs = acs'; end
    %bnd = (sqrt(2)*erfinv(0.95))./sqrt(T);
    %idx = find(abs(acs)>bnd);
    %isit       = abs(acs)>bnd & (1:T);
    %where2stop = find(isit==0); %finds the break point -- intercept 
    %BE CAREFUL:
    %this here is different from Toeplitz ME version, because here we don't
    %have the 0lag, but in that setting the 0lag is there and it is the
    %diag. 
    where2stop = FindBreakPoint(acs,T);
    
    if ~where2stop %if there was nothing above the CI...
        srnkd_ts = zeros(1,T);
    else
        % where2stop = where2stop(1)-1; %-1 because we want to stop before intercept
        % srnkd_ts   = tukeytaperme(ts,where2stop);
        srnkd_ts   = curbtaperme(acs,T,where2stop);
    end
end
%--------------------------------------------------------------------------
function where2stop = FindBreakPoint(acs,T)
% this finds the breaking points for shrinking the AC. 
% Nothing serious, just might help with speed...
% SA, Ox, 2018
    if ~sum(ismember(size(acs),T)); error('There is something wrong, mate!'); end
    if size(acs,2) ~= T; acs = acs'; end
    
    bnd        = (sqrt(2)*erfinv(0.95))./sqrt(T);
    %idx        = find(abs(acs)>bnd);
    isit       = abs(acs)>bnd & (1:T);
    where2stop = find(isit==0); %finds the break point -- intercept 
    
    if where2stop(1)==1
        where2stop = 0; 
    else
        where2stop = where2stop(1)-1; 
    end;
end
%--------------------------------------------------------------------------
function ct_ts=curbtaperme(acs,T,M)
% Curb the autocorrelations, according to Anderson 1984.
% Multi-dimensional, and therefore is fine!
%SA, Ox, 2018
    if ~sum(ismember(size(acs),T)); error('There is something wrong, mate!'); end
    if size(acs,2) ~= T; acs = acs'; end
    
    M          = round(M);
    msk        = zeros(size(acs));
    msk(:,1:M) = 1;
    ct_ts      = msk.*acs;
end
%--------------------------------------------------------------------------
function tt_ts=tukeytaperme(acs,T,M)
%performs Single Tukey Tapering for given length of window, M, and initial
%value, intv. intv should only be used on crosscorrelation matrices.
%
%NB! There used to be initialisation parameters here before, intv. I 
%remoeved it because we now start with the second elements of the ACF anyways. 
%
%SA, Ox, 2018
    if ~sum(ismember(size(acs),T)); error('There is something wrong, mate!'); end
    if size(acs,2) ~= T; acs = acs'; end
    %if ~exist('intv','var'); intv = 1; warning('Oi!'); end;
    M          = round(M);
    tt_ts      = zeros(size(acs));
    tt_ts(1)   = 1;
    tt_ts(2:M) = (1 + cos([2:M] .* pi ./ M))./2 .* acs(2:M);
    %figure; plot(tt_ts); hold on; plot(acs); 
end
