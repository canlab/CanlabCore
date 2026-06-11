function [uncorr_np_pval,uncorr_param_pval] = NonParamCorrInf(ts,T,nRlz,Taper,TaperMethod,TaperParam)

%--------------------------------------------------------------
% the input section looks fairly shit, change it later to more user
% friendly settings:

if size(ts,1) ~= 2 || size(ts,2) ~= T
    error('PFF!')
end

if ~exist('Taper','var') || ~strcmpi(Taper,'taper')
    Taper = ''; 
end;

if ~exist('TaperMethod','var') || ~any(strcmpi(TaperMethod,{'shrink','tukey','curb'})) 
    TaperMethod = ''; 
end;

if ~exist('TaperParam','var')
    TaperParam = []; 
end;

if ~exist('nRlz','var') || isempty(nRlz)
    nRlz = 5000;
end

disp([num2str(nRlz) '-' Taper ',' TaperMethod ',' num2str(TaperParam)])
%--------------------------------------------------------------
[~,stat0]=MonsterEquation(ts,T,Taper,TaperMethod,TaperParam);
Z_fnp = stat0.z.rzf;
%Z_fisher_np = atanh(corr(ts(1,:)',ts(2,:)')).*sqrt(T-3);
%-----
acx = AC_fft(ts(1,:),T);
acy = AC_fft(ts(2,:),T);

AC_X = MakeMeCovMat(acx(2:end),T); Kx = sqrtm(AC_X);
AC_Y = MakeMeCovMat(acy(2:end),T); Ky = sqrtm(AC_Y);

z_ME_np     = zeros(1,nRlz); 
%z_fisher_np = zeros(1,nRlz);
for i=1:nRlz    
    %-This is identical to corrautocorr([0 0],0,0,T), but much faster!
    Wx = randn(1,T)'; 
    Wy = randn(1,T)';
    Gx = Kx*Wx;     
    Gy = Ky*Wy;
    %-------------------------------------------------
    %-Remove this after the paper, this eats half of the exec time! 
    %z_fisher_np(i) = atanh(corr(Gx,Gy)).*sqrt(T-3);
    %---- 
    [~,stat] = MonsterEquation([Gx,Gy]',T,Taper,TaperMethod,TaperParam);
    z_ME_np(i) = stat.z.rzf;
    clear stat Gx Gy Wx Wy
end
uncorr_np_pval = sum(z_ME_np>Z_fnp)./nRlz;

%-outputs----
%uncorr_np_pval_fisher = sum(z_fisher_np>Z_fisher_np)./nRlz;
if ~uncorr_np_pval % we don't have zero pvalues, realisations weren't enough dummy!
    uncorr_np_pval = eps;
end
uncorr_param_pval = stat0.p.f_Pval;
