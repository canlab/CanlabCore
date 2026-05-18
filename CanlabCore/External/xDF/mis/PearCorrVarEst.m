function [ASAt,Stat]=PearCorrVarEst(ts,T,varargin)
% Estimates the Monster Equation!
%
%%%INPUTS:
%   ts: Time series as a 2D matrix. 
%   T : Number of data-points. Just for the sake of sanity checks
%
%   Optionals:
%   'taper' : uses a tapering method to denoise AC functions
%   'TVOff' : if an estimate exceeed the theoritical variance of a white
%   noise then it curbs the estimate back to (1-rho^2)^2/T. If you want it
%   off, trigger 'TVOff'
%%%OUTPUTS:
%   ASAt : a 2D matrix of size IxI with diagonal set to zero
%   Stat : Contains every other non very important stuff!
%%%DEPENDECIES:
%   AC_fft.m : estimates ACF super quick via FFT
%   xC_fft.m : estimates cross corr functions super quick via FFT
%
%%%REFERENCES:
%   Variance of Pearson's correlations of serially-correlated time series
%   Soroosh Afyouni & Thomas E. Nichols
%   2018
%   University of Oxford
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2018
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________
    
%-----------------------------------------------------------
error('DO NOT USE ME, USE xDF.m')
%-----------------------------------------------------------

    if  size(ts,2) ~= T
        ts = ts';
        warning('Oi!')
    end
    
    W2S = []; TVflag = 1;
    
    nn  = size(ts,1);
    %ts  = dtrend(ts);
    ts  = ts./std(ts,[],2); %standardise
    %Corr----------------------------------------------------------------------
    rho   = corr(ts');
    rho(1:nn+1:end) = 0;
    %Autocorr------------------------------------------------------------------
    [ac] = AC_fft(ts,T); %demean the time series
    ac   = ac(:,2:T-1);
    
    nLg  = T-2;       %if lag0 was the 0th element. Also, the ACF has T-1 dof. eh?  
    %Cross-corr---------------------------------------------------------------- 
    xcf = xC_fft(ts,T);
    acx_n      = flip(xcf(:,:,2:T-1),3);
    acx_p      = xcf(:,:,T+1:end-1);
    
    %figure; plot(ac')
    %figure; plot(squeeze(acx_n(1,2,:))); hold on; plot(squeeze(acx_p(1,2,:)));  
    %----MEMORY SAVE----
    clear ts 
    %-------------------
        
    
    if sum(strcmpi(varargin,'TVOff'))
        %disp('Anything beyond theoritical variance is forced down to theoritical variance (1-rho^2)^2.')
        TVflag = 0;
    end    
    
    if sum(strcmpi(varargin,'taper'))
        mth = varargin{find(strcmpi(varargin,'taper'))+1};
        if strcmpi(mth,'tukey')
    %Tukey Tappering----------------------------------------
            M = round(varargin{find(strcmpi(varargin,'tukey'))+1});
            if isempty(M); error('you MUST set a tukey factor.'); end;

            %disp([mth 'ed with ' num2str(M) ' length was used.'])
            for in=1:nn    
                ac(in,:) = tukeytaperme(ac(in,:),nLg,M);
                for jn=1:nn 
                    acx_n(in,jn,:) = tukeytaperme(squeeze(acx_n(in,jn,:)),nLg,M);
                    acx_p(in,jn,:) = tukeytaperme(squeeze(acx_p(in,jn,:)),nLg,M);
                end
            end
            
    %Shrinking------------------------------------------------
        elseif strcmpi(mth,'shrink')
            M = round(varargin{find(strcmpi(varargin,'shrink'))+1});
            
            %disp('TAPER SHRINK')
            
            if M>1
                error('Not yet, mate!')
%                 for in=1:nn    
%                     ac(in,:)  = ShrinkPeriod(ac(in,:),M);
%                     for jn=1:nn 
%                         acx_n(in,jn,:) = ShrinkPeriod(acx_n(in,jn,:),M);
%                         acx_p(in,jn,:) = ShrinkPeriod(acx_p(in,jn,:),M);
%                     end
%                 end
            elseif M==1
                
                
                for in=1:nn
                    for jn=1:nn
                        W2S(in,jn) = max([FindBreakPoint(ac(in,:),nLg) FindBreakPoint(ac(jn,:),nLg)]);
                    end
                end
                
                for in=1:nn    
                    ac(in,:)= shrinkme(ac(in,:),nLg);
%                     for jn=1:nn 
%                         acx_n(in,jn,:) = shrinkme(squeeze(acx_n(in,jn,:)),nLg);
%                         acx_p(in,jn,:) = shrinkme(squeeze(acx_p(in,jn,:)),nLg);
%                     end
                    for jn=1:nn 
                        acx_n(in,jn,:) = curbtaperme(squeeze(acx_n(in,jn,:)),nLg,W2S(in,jn));
                        acx_p(in,jn,:) = curbtaperme(squeeze(acx_p(in,jn,:)),nLg,W2S(in,jn));
                    end

                end  
            else
                error('What are you up to mate?!')
            end
            %figure; plot(squeeze(acx_n(1,2,:))); hold on; plot(squeeze(acx_p(1,2,:)));             
            %figure; plot(ac(1,:)); hold on; plot(ac(2,:));
    %Curbing------------------------------------------------
        elseif strcmpi(mth,'curb')
            M = round(varargin{find(strcmpi(varargin,'curb'))+1});

            for in=1:nn    
                ac(in,:) = curbtaperme(ac(in,:),nLg,M);
                for jn=1:nn 
                    %size(squeeze(acx_n(in,jn,:))), nLg
                    acx_n(in,jn,:) = curbtaperme(squeeze(acx_n(in,jn,:))',nLg,M);
                    acx_p(in,jn,:) = curbtaperme(squeeze(acx_p(in,jn,:))',nLg,M);
                end
            end
    %--------------------------------------------------------------------------
        else
            error('choose one of these; shrink | tukey | curb as tapering option.')
        end
    end
   
    
    
%Crazy, eh?
wgt     = (nLg:-1:1);
wgtm3   = reshape(repmat((repmat(wgt,[nn,1])),[nn,1]),[nn,nn,numel(wgt)]); %this is shit, eats all the memory!
wgtm2   = repmat(wgt,[nn,1]);
Tp      = T-1;

ASAt = [Tp                  .* (1-rho.^2).^2 ...
       + rho.^2             .* sum(SumMat((wgtm2.*ac.^2),nLg),3) ...                % 1    -- AC 
       + 2 .* wgt           .* ac*ac'...                                            % 5    -- AC
       + 2 .* (rho.^2 + 1)  .* sum((wgtm3.*acx_n.*acx_p),3) ...                     %2 & 3 -- XC
       - 2 .* rho           .* sum(wgtm3.*SumMat(ac,nLg).*(acx_n+acx_p),3)]./(T.^2);  %4 -- This this the only term which we can't seperate the AC and XC!

%----MEMORY SAVE----
clear wgtm3 acx_* ac 
%-------------------

%Keep your wit about you!
TV = (1-rho.^2).^2./T;
%TV
if sum(sum(ASAt < TV)) && TVflag
    % Considering that the variance can *only* get larger in presence of autocorrelation.   
    idx_ex       = find(ASAt < TV);
    ASAt(idx_ex) = TV(idx_ex);
    
    warning([num2str(numel(idx_ex)) ' edges had variance smaller than the textbook variance!'])
end  
Stat.TV    = TV;

% diagonal is rubbish;
ASAt(1:nn+1:end) = 0;

%------- Test Stat-----------------------
%Pearson's turf -- We don't really wanna go there, eh?
%rz      = rho./sqrt((ASAt));     %abs(ASAt), because it is possible to get negative ASAt!
%r_pval  = 2 * normcdf(-abs(rz)); %both tails
%r_pval(1:nn+1:end) = 0;          %NaN screws up everything, so get rid of the diag, but becareful here. 

%Our turf--------------------------------
rf      = atanh(rho);
sf      = ASAt./((1-rho.^2).^2);    %delta method; make sure the N is correct! So they cancelled out?!
rzf     = rf./sqrt(sf);
rzf(1:nn+1:end) = 0;
f_pval  = 2 .* normcdf(-abs(rzf));  %both tails
f_pval(1:nn+1:end) = 0;             %NaN screws up everything, so get rid of the diag, but becareful here. 

Stat.z.rzf = rzf;
Stat.p.f_Pval = f_pval;

%Fisher's turf---------------------------
% rf          = atanh(rho);
% edf         = 1./ASAt;                          %Effective Degrees of Freedom
% rzfish      = rf.*sqrt(edf-3);
% rzfish(1:nn+1:end)  = 0;
% f_pval_fish         = 2 .* normcdf(-abs(rzfish));  %both tails
% f_pval_fish(1:nn+1:end) = 0;                %NaN screws up everything, so get rid of the diag, but becareful here. 
% 
% Stat.z.rzfish = rzfish;
% Stat.p.f_PvalFish = f_pval_fish;
%-------Stat-----------------------------
Stat.W2S = W2S;
end
%--------------------------------------------------------------------------          

%--------------------------------------------------------------------------
function [SM0] = SumMat(Y0,T)
    %hopefully faster than Matlab's dumb sum of sum!
    %SA, Ox, 2018

    %becareful with this function, this is really tricky to use!
    if ~sum(ismember(size(Y0),T)); error('There is something wrong, mate!'); end
    if size(Y0,1) ~= T; Y0 = Y0'; end

    nn  = size(Y0,2);
    idx = find(triu(ones(nn),1))';
    %SM  = zeros(nn);
    SM0 = zeros(nn,nn,T);
    for i=idx
        [x,y]      = ind2sub(nn,i);
        SM0(x,y,:) = (Y0(:,x)+Y0(:,y));
        SM0(y,x,:) = (Y0(:,y)+Y0(:,x));
    end
    %SM = sum(SM0,3);
end

%--------------------------------------------------------------------------

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
% Curb the autocorrelations, according to Anderson 1984
% multi-dimensional, and therefore is fine!
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
    %tt_ts(:,1) = intv;
    tt_ts(1:M) = (1+cos([1:M].*pi./M))./2.*acs(1:M);
    %figure; plot(tt_ts); hold on; plot(acs); 
end

% function sY=ShrinkPeriod(Y,WhichPeak)
% %WhichPeak indetified up until how many bumps we should continue (above-CI
% %bumps, of course).
% %NB! This function uses pieces from an external package.
%
% %SA, Ox, 2018
%
%     T = numel(Y);
%     bnd=(sqrt(2)*erfinv(0.95))./sqrt(T);
%     
%     P  = round(InterX([1:T;Y],[1:T;zeros(1,T)]));
%     P  = P(1,:); 
%     P0 = [1 P(1,:)];
%     
%     if WhichPeak>numel(P)
%         warning('There are less peaks than you asked for.')
%         WhichPeak=numel(P);
%     end
%     
%     Idx=zeros(1,WhichPeak);
%     for p=1:WhichPeak; pY = Y(P0(p):P(p)); if any(abs(pY)>bnd); Idx(p)=1; end; end;
%     
%     for i=Idx
%         if ~i; break; end; 
%         Idx(i)=1; 
%     end
%     
%     sY0 = Y(1:max(P(find(Idx))));
%     sY  = zeros(1,T);
%     sY(1:numel(sY0)) = sY0;
% end
% function P = InterX(L1,L2)
%        
%     %...Preliminary stuff
%     x1  = L1(1,:)';  x2 = L2(1,:);
%     y1  = L1(2,:)';  y2 = L2(2,:);
%     dx1 = diff(x1); dy1 = diff(y1);
%     dx2 = diff(x2); dy2 = diff(y2);
%     
%     %...Determine 'signed distances'   
%     S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
%     S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
%     
%     C1 = feval(@le,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
%     C2 = feval(@le,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';
% 
%     %...Obtain the segments where an intersection is expected
%     [i,j] = find(C1 & C2); 
%     if isempty(i),P = zeros(2,0);return; end;
%     
%     %...Transpose and prepare for output
%     i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
%     L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
%     i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
%     
%     %...Solve system of eqs to get the common points
%     P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
%                 dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
%               
%     function u = D(x,y)
%         u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
%     end
% end

