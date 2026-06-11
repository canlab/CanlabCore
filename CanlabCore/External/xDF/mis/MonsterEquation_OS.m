function [ASAt,CnR] = MonsterEquation_OS(rho,ACx,ACy,T)
% [ASAt,CnR] = MonsterEquation_OS(rho,ACx,ACy,T)
% CALCULATE the variance for given rho ACx ACy and T
% ONLY use this for oracle simulations. 
%
%%%INPUTS: all first three inputs should be the *true* values. 
%
%%%NOTES:
% 
%   See Algorithm 1 of SuppMat in the paper. 
%
%%%REFERENCES:
%  Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. 
% "Effective Degrees of Freedom of the Pearson's Correlation Coefficient 
%  under Serial Correlation." bioRxiv (2018): 453795.
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2018
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

% M = 1; 
% if M>1
%     ACx  = ShrinkPeriod(ACx,M);
%     ACy  = ShrinkPeriod(ACy,M);
%     elseif M==1
%     ACx  = shrinkme(ACx);
%     ACy  = shrinkme(ACy);  
% end

Sigma_x  = MakeMeCovMat(ACx,T,0);
Sigma_y  = MakeMeCovMat(ACy,T,0);

Kx = sqrtm(Sigma_x);
Ky = sqrtm(Sigma_y);

d   = diag(Kx*Ky');

as  = min(d);
rho = rho.*as;

Sigma_xy = rho.*((Kx*Ky')./as);

%figure; imagesc(Sigma_xy)

%Sigma_xy(1:20,1:20)

% figure; plot(d); hold on; plot(diag((Kx*Ky')./as)); hold off; 
% figure; 
% subplot(3,1,1); hold on; axis tight; 
% imagesc(Sigma_x)
% subplot(3,1,2); hold on; axis tight; 
% imagesc(Sigma_y)
% subplot(3,1,3); hold on; axis tight; 
% imagesc(Sigma_xy)

%Monster Eq -----------------------------------------------
% SigX2    =  trace(Sigma_x ^2);
% SigY2    =  trace(Sigma_y ^2);
% SigXSigY =  trace(Sigma_x * Sigma_y);
% SigXY2    = trace(Sigma_xy^2);
% SigXSigXY = trace(Sigma_x * Sigma_xy);
% SigYSigXY = trace(Sigma_y * Sigma_xy);
% ASAt=((rho.^2./2).* SigX2... 
%     + (rho.^2./2) .* SigY2...
%     + rho.^2      .* SigXY2...
%     - 2.*rho      .* SigXSigXY...
%     - 2.*rho      .* SigYSigXY... 
%     + SigXSigY...
%     + SigXY2)/T.^2;

%----------------------------------------------------------
%xDF
ASAt = ((rho.^2./2) .* (trace(Sigma_x^2) + trace(Sigma_y^2))...                     % 1 & 2
        + trace(Sigma_x*Sigma_y) + rho^2 .* trace(Sigma_xy*Sigma_xy')...            % 3&4
        - rho .* trace(Sigma_x*Sigma_xy) - rho .* trace(Sigma_x*Sigma_xy')...       % 5&6
        - rho .* trace(Sigma_y*Sigma_xy) - rho .* trace(Sigma_y*Sigma_xy')...       % 7&8
        + trace(Sigma_xy^2))./T^2; %9

%----------------------------------------------------------
SigXSigY =  trace(Sigma_x * Sigma_y);
CnR = SigXSigY/T^2;
%(1-rho.^2).^2/T

end

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
% function ct_ts=curbtaperme(ts,M)
% % Curb the autocorrelations, according to Anderson 1984
%     M          = round(M);
%     msk        = zeros(size(ts));
%     msk(:,1:M) = 1;
%     ct_ts      = msk.*ts;
% end
% function tt_ts=tukeytaperme(ts,M,intv)
% %performs Single Tukey Tapering for given length of window, M, and initial
% %value, intv. intv should only be used on crosscorrelation matrices. 
%     if ~exist('intv','var'); intv = 1; warning('Oi!'); end;
%     M          = round(M);
%     tt_ts      = zeros(size(ts));
%     tt_ts(:,1) = intv;
%     tt_ts(2:M) = (1+cos([2:M].*pi./M))./2.*ts(2:M);
% end

% function sY=ShrinkPeriod(Y,WhichPeak)
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



