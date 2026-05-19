function Wy = WhiteningCorr(Y,L,varargin)
% Wy = Whitening(Y,L,varargin)
% Decorelate and Whiten the data by:
% PCA:      Wy = WhiteningCorr(Y,L,'Method','PCA')
% ZCA:      Wy = WhiteningCorr(Y,L,'Method','ZCA')
% Cholesky: Wy = WhiteningCorr(Y,L,'Method','Chol')
% SQRTM:    Wy = WhiteningCorr(Y,L,'Method','sqrtm')
%
%   NB! Becareful about number of time series, the estm gets crappier as 
%   the number of time series grows.  
%
%%%REFERENCES:
%  Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. 
% "Effective Degrees of Freedom of the Pearson's Correlation Coefficient 
%  under Serial Correlation." bioRxiv (2018): 453795.
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________
%
if size(Y,2)~=L
    Y=Y'; %IxT
end
% I = size(Y,1);

dmY   = Y-mean(Y,2);
sigma = dmY*dmY'./(L-1);


if nargin>2 && ~sum(strcmpi(varargin,{'Method','Decorrelate'}))
    error('Choose either Method or Decorrelate.')
end

if sum(strcmpi(varargin,'Method'))
    MethodStr = varargin{find(strcmpi(varargin,'Method'))+1};
    if strcmpi(MethodStr,'PCA') 
        %disp('**PCA**')
        [U,S,~] = svd(sigma);                        
        Wy = diag(1./sqrt(diag(S) + eps)) * U' * dmY;
    elseif strcmpi(MethodStr,'ZCA')
        %disp('**ZCA**')
        [U,S,~] = svd(sigma);
        Wy = U * diag(1./sqrt(diag(S) + eps)) * U' * dmY;
    elseif strcmpi(MethodStr,'Chol') || strcmpi(MethodStr,'Cholesky')
        R  = chol(sigma)';
        Wy = pinv(R)*dmY; %FUCK!
    elseif strcmpi(MethodStr,'SQRTM')
        %disp('**SQRTM**')
        R  = sqrtm(sigma)';
        Wy = pinv(R)*dmY; %FUCK!         
    end
elseif sum(strcmpi(varargin,'Decorrelate')) % this only decorrelates!
    [U,~,~] = eig(sigma);    
    Wy = (dmY'*U')';
else
    R  = sqrtm(sigma)';
    Wy = pinv(R)*dmY; 
   
    
end

%corr(PCA')
%corr(ZCA')
% corr(Wy')
% 
% dmY_ac1=autocorr(dmY(1,:),L-1);
% Wy_ac1=autocorr(Wy(1,:),L-1);
% 
% dmY_ac2=autocorr(dmY(2,:),L-1);
% Wy_ac2=autocorr(Wy(2,:),L-1);
% 
% figure; hold on; box on;
% subplot(2,1,1); hold on;
% plot(dmY_ac1,'linewidth',4,'color',[.5 .5 .5 .5])
% plot(Wy_ac1,'linewidth',1.1,'color','r')
% subplot(2,1,2); hold on;
% plot(dmY_ac2,'linewidth',4,'color',[.5 .5 .5 .5])
% plot(Wy_ac2 ,'linewidth',1.1,'color','r')

% legend({'orig','Prewhitened'})
% figure; hold on; box on;
% scatter(dmY(1,:),dmY(2,:),'r.')
% scatter(Wy(1,:),Wy(2,:),'k.')

