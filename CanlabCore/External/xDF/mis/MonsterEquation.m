function [ASAt,Stat]=MonsterEquation(ts,T,varargin)

    if size(ts,1) ~= 2 || size(ts,2) ~= T
        error('PFF!')
    end
    %ts  = dtrend(ts);
    ts  = ts./std(ts,[],2); %standardise
    %Corr----------------------------------------------------------------------
    c   = corr(ts');
    rho = (c(1,2));
    %Autocorr------------------------------------------------------------------
    [ac] = AC_fft(ts,T); %demean the time series
    ac_x = ac(1,1:T-1);
    ac_y = ac(2,1:T-1);
    %Cross-corr---------------------------------------------------------------- 
    [xcf,lags]  = crosscorr(ts(1,:),ts(2,:),T-1); %demean the time series
    acx_n = fliplr(xcf(2:T));
    acx_p = xcf(T:end-1);

    
%     figure; hold on;  
%     plot(acx_n); plot(acx_p)
    
    if sum(strcmpi(varargin,'taper'))
        mth = varargin{find(strcmpi(varargin,'taper'))+1};
        if strcmpi(mth,'tukey')
    %Tukey Tappering----------------------------------------
            M = round(varargin{find(strcmpi(varargin,'tukey'))+1});
            if isempty(M); error('you MUST set a tukey factor.'); end;

            %disp([mth 'ed with ' num2str(M) ' length was used.'])

            ac_x = tukeytaperme(ac_x,M);
            ac_y = tukeytaperme(ac_y,M);

            %acx_n = shrinkme(acx_n);
            %acx_p = shrinkme(acx_p);

            acx_n = tukeytaperme(acx_n,M);
            acx_p = tukeytaperme(acx_p,M);
    %Shrinkage------------------------------------------------
        elseif strcmpi(mth,'truncate')
            M = round(varargin{find(strcmpi(varargin,'truncate'))+1});
            
            %disp('TAPER SHRINK')
            
            if M>1
                ac_x  = ShrinkPeriod(ac_x,M);
                ac_y  = ShrinkPeriod(ac_y,M);

                acx_n = ShrinkPeriod(acx_n,1);
                acx_p = ShrinkPeriod(acx_p,1);
                
            elseif M==1
                [ac_x,w2s0] = shrinkme(ac_x,T-1);
                [ac_y,w2s1] = shrinkme(ac_y,T-1);
                
                w2s = max([w2s0 w2s1]);
                
                %figure; plot(ac_x); hold on; plot(ac_y); 
                
                acx_n = curbtaperme(acx_n,w2s);
                acx_p = curbtaperme(acx_p,w2s);
                %acx_n = shrinkme(acx_n,T-1);
                %acx_p = shrinkme(acx_p,T-1);    
                
                %figure; plot(acx_n); hold on; plot(acx_p); 
            else
                error('What are you up to mate?!')
            end
    %Curbing------------------------------------------------
        elseif strcmpi(mth,'curb')
            M = round(varargin{find(strcmpi(varargin,'curb'))+1});

            ac_x = curbtaperme(ac_x,M);
            ac_y = curbtaperme(ac_y,M);

            acx_n = curbtaperme(acx_n,M);
            acx_p = curbtaperme(acx_p,M);
    %--------------------------------------------------------------------------
        else
            error('choose shrink, tukey and curb as taper option.')
        end
    end
    %plot(acx_n)
    
    Sigma_x  = toeplitz(ac_x);
    Sigma_y  = toeplitz(ac_y);
    Sigma_xy = triu(toeplitz(acx_n))+tril(toeplitz(acx_p),-1);

    %-------ME
    SigX2     = trace(Sigma_x ^2);
    SigY2     = trace(Sigma_y ^2);
    SigXSigY  = trace(Sigma_x * Sigma_y);
    SigXY2    = trace(Sigma_xy^2);
    SigXSigXY = trace(Sigma_x * Sigma_xy);
    SigYSigXY = trace(Sigma_y * Sigma_xy);

    ASAt      = ((rho.^2./2) .* SigX2... 
                +(rho.^2./2) .* SigY2...
                +rho.^2      .* SigXY2...
                -2.*rho      .* SigXSigXY...
                -2.*rho      .* SigYSigXY... 
                +SigXSigY...
                +SigXY2)./T.^2;
            
    Stat.ME.trSigX2       = SigX2;
    Stat.ME.trSigY2       = SigY2;
    Stat.ME.trSigXSigY    = SigXSigY;
    Stat.ME.trSigXY2      = SigXY2; 
    Stat.ME.trSigXSigXY   = SigXSigXY;
    Stat.ME.trSigYSigXY   = SigYSigXY; 
    Stat.CnR=SigXSigY./T^2; %just to check how bad the others are doing!


%%%%%%%THIS THIS BULLSHT!    
%ASAt    
%-----ME, but damn faster this time!
% nLg     = T-2;        %if lag0 was the 0th element. Also, the ACF has T-1 dof. eh?  
% 
% wgt     = (nLg:-1:1);
% Tp      = T-1;
% 
% LAMBDAX = ac_x (2:end);
% LAMBDAY = ac_y (2:end);
% RHOp    = acx_p(2:end);
% RHOn    = acx_n(2:end);

% size(LAMBDAX),size(LAMBDAY),size(RHOp),size(RHOn)

% ASAt = [Tp*(1-rho.^2).^2 ...
%     + rho.^2 * sum( wgt .* (LAMBDAX.^2 + LAMBDAY.^2 + 2*RHOp.*RHOn))...
%     - 2 * rho* sum( wgt .* (RHOp+RHOn) .* (LAMBDAX+LAMBDAY))...
%     + 2 *      sum( wgt .* (RHOp.*RHOn  + LAMBDAX.*LAMBDAY))]./T.^2;
%%%%%%%THIS THIS BULLSHT!  


%Keep your wit about you!
TV = (1-rho.^2).^2./T;

%ASAt

% if ASAt<TV
%     ASAt=TV; 
% end; 
    %------- Test Stat
    %Pearson's turf
    rz     = rho./sqrt((ASAt));    %abs(ASAt), because it is possible to get negative ASAt!
    r_pval = 2 * normcdf(-abs(rz));  %both tails
    %Fisher's turf
    rf   = atanh(rho);
    sf   = ASAt./((1-rho^2)^2); %delta method
    rzf  = rf./sqrt(sf);
    f_pval = 2 * normcdf(-abs(rzf)); %both tails
    %-------Stat

    Stat.p.r_Pval = r_pval;
    Stat.p.f_Pval = f_pval;

    Stat.z.rz  = rz;
    Stat.z.rzf = rzf;
end
%--------------------------------------------------------------------------          

% function srnkd_ts=shrinkme(ts)
% %Shrinks the *early* bucnhes of autocorr coefficients beyond the CI.
%     L = numel(ts);
%     bnd = (sqrt(2)*erfinv(0.95))./sqrt(L);
%     idx = find(abs(ts)>bnd);
%     isit       = abs(ts)>bnd & (1:L);
%     where2stop = find(isit==0);
%     where2stop = where2stop(1);
%     % srnkd_ts   = tukeytaperme(ts,where2stop);
%     srnkd_ts   = curbtaperme(ts,where2stop);
% end
function [srnkd_ts where2stop]=shrinkme(ts,T)
%Shrinks the *early* bucnhes of autocorr coefficients beyond the CI.
    if ~sum(ismember(size(ts),T)); error('There is something wrong, mate!'); end
    if size(ts,2) ~= T; ts = ts'; end
        
    bnd = (sqrt(2)*erfinv(0.95))./sqrt(T);
    
    idx = find(abs(ts)>bnd);
    isit       = abs(ts)>bnd & (1:T);
    where2stop = find(isit==0);
    where2stop = where2stop(1)-1; %-1 because we want to stop before intercept
    % srnkd_ts   = tukeytaperme(ts,where2stop);
    srnkd_ts   = curbtaperme(ts,where2stop);
end
function ct_ts=curbtaperme(ts,M)
% Curb the autocorrelations, according to Anderson 1984
    M          = round(M);
    msk        = zeros(size(ts));
    msk(:,1:M) = 1;
    ct_ts      = msk.*ts;
end
function tt_ts=tukeytaperme(ts,M)
%performs Single Tukey Tapering for given length of window, M, and initial
%value, intv. intv should only be used on crosscorrelation matrices. 
%
%NB! There used to be initialisation parameters here before, intv. I 
%remoeved it because we now start with the second elements of the ACF anyways. 

    %if ~exist('intv','var'); intv = 1; warning('Oi!'); end;
    M          = round(M);
    tt_ts      = zeros(size(ts));
    %tt_ts(:,1) = intv;
    tt_ts(1:M) = (1+cos([1:M].*pi./M))./2.*ts(1:M);
end

function sY=ShrinkPeriod(Y,WhichPeak)

    T = numel(Y);
    bnd=(sqrt(2)*erfinv(0.95))./sqrt(T);
    
    P  = round(InterX([1:T;Y],[1:T;zeros(1,T)]));
    P  = P(1,:); 
    P0 = [1 P(1,:)];
    
    if WhichPeak>numel(P)
        warning('There are less peaks than you asked for.')
        WhichPeak=numel(P);
    end
    
    Idx=zeros(1,WhichPeak);
    for p=1:WhichPeak; pY = Y(P0(p):P(p)); if any(abs(pY)>bnd); Idx(p)=1; end; end;
    
    for i=Idx
        if ~i; break; end; 
        Idx(i)=1; 
    end
    
    sY0 = Y(1:max(P(find(Idx))));
    sY  = zeros(1,T);
    sY(1:numel(sY0)) = sY0;
end
function P = InterX(L1,L2)
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(@le,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(@le,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end

% figure; 
% subplot(3,2,1); hold on;
% title('Sigma_{xy}')
% imagesc(Sigma_xy,[-1 1])
% disp(['Sigma_{xy}: ' num2str(trace(Sigma_xy))])
% 
% subplot(3,2,2); hold on; 
% title('Sigma_x x Sigma_y')
% imagesc(Sigma_x * Sigma_y,[-1 1])
% disp(['Sigma_{x}Sigma_{y}: ' num2str(trace(Sigma_x * Sigma_y)) ]);
% 
% subplot(3,2,3); hold on; 
% title('Sigma_x x Sigma_{xy}')
% imagesc(Sigma_x * Sigma_xy,[-1 1])
% disp(['Sigma_{x}Sigma_{xy}: ' num2str(trace(Sigma_x * Sigma_xy)) ]);
% 
% subplot(3,2,4); hold on; 
% title('Sigma_y x Sigma_{xy}')
% imagesc(Sigma_y * Sigma_xy,[-1 1])
% disp(['Sigma_{y}Sigma_{xy}: ' num2str(trace(Sigma_y * Sigma_xy))])
% 
% subplot(3,2,5); hold on; 
% title('Sigma x')
% imagesc(Sigma_x,[-1 1])
% disp(['Sigma_{x}: ' num2str(trace(Sigma_x))])
% 
% subplot(3,2,6); hold on; 
% title('Sigma y')
% imagesc(Sigma_y,[-1 1])
% disp(['Sigma_{y}: ' num2str(trace(Sigma_y))])

% CnRe = trace(Sigma_x * Sigma_y); % Clifford and Richardson; under the null
% 
% ASA=((rho.^2./2).* trace(Sigma_x ^2        )...
%     +rho.^2     .* trace(Sigma_xy^2        )...
%     -2.*rho     .* trace(Sigma_x * Sigma_xy)...
%     +(rho.^2./2).* trace(Sigma_y ^2        )...
%     -2.*rho     .* trace(Sigma_y * Sigma_xy)...
%     + CnRe...
%     + trace(Sigma_xy^2))./T^2;
