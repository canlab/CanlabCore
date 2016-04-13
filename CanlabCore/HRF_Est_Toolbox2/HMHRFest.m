function [Res] = HMHRFest(y, Runs, TR, nbasis, norder)
% HRF estimation algorithm
%
% :Inputs:
%
%   **y:**
%        Data matrix (#time points) by (#subjects) by (#voxels)
%
%   **Runs:**
%        Stick functions for each subject (#time points) by (#conditions) by (#subjects) 
%
%   **TR:**
%        Time resolution
%
%   **nbasis:**
%        Number of b-spline basis
%
%   **norder:**
%        Order of b-spline basis
%
% ..
%    By David Degras and Martin Lindquist
%    Created: 08/22/12
%    Last edited: 03/27/14
% ..


[len, sub, voxpar] = size(y);
[~, L, ~] = size(Runs);

% Set initial values
%
% TR = 1;                 % TR of experiment
% nbasis = 20;            % Number of b-spline basis
% norder = 6;             % Order of b-spline basis

I = 1;                  % Number of groups
N1 = sub;               % N1: Number of subjects in Group 1
Tlen = 30/TR;           % Length of HRF
q = 2;                  % Number of nuisance parameters
p = 1;                  % AR order
lambda = 50;            % Smoothing parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Info about the parcellation

voxels = sum(voxpar);               % Total number of voxels 
parcels = length(voxpar);           % Number of parcels
parind = zeros(voxpar(1),1)+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a struct containing all initial values

Init = {};
Init.I = I;
Init.sub = sub;
Init.len = len;
Init.N1 = N1;
Init.L = L;
Init.TR = TR;
Init.Tlen = Tlen;
Init.nbasis = nbasis;
Init.norder = norder;
Init.q = q;
Init.p = p;
Init.lambda = lambda;
Init.voxels = voxels;
Init.parcels = parcels;
Init.voxpar = voxpar;
Init.parind = parind;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate basis sets for potential HRF shapes.

hphi = zeros(3,nbasis);        % Basis sets for 3 potenital HRF shapes
hf = CanonicalBasisSet(TR);

% Create bspline basis set
basis = create_bspline_basis([0,Tlen], nbasis+8, norder);    
B = eval_basis((1:Tlen),basis);
B = B(:,6:end-3);


Init.B = B;
SM = zeros(Tlen,3);
for i=1:3,

    Run = zeros(Tlen,1);
    Run(1:(2*(i-1)+1)) = 1;
    s = conv(hf,Run);
    s = s(1:Tlen);
    s = s./max(s);

    SM(:,i) = s;
    tmp = (inv(B'*B)*B'*s)';
    hphi(i,:) = tmp./sum(tmp.^2);

end

hphi = orth(hphi')';            % Orthogonalize

Init.hphi = hphi;


G0 = (inv(B'*B)*B'*hf);         % Need these values for Hotellings test.



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of HRF and nuisance parameters

Beta = zeros(L, voxels);                    % Voxel-wise beta
G = zeros(nbasis, voxels);                  % Voxel-wise gamma
K = zeros(N1*nbasis*L, voxels);             % Voxel-wise ksi
% E = zeros(len, voxels);                   % Voxel-wise epsilon
Ar = zeros(p,voxels);                       % Voxel-wise AR coeficients
S2e = zeros(1,voxels);                      % Voxel-wise withing-subject variance
CC = zeros(1,voxels);                       % Voxel-wise Lagrange Multiplier

% Create initial covariance matrix (indentity matrix for all subjects)

iV = zeros(len,len,sub);
for j=1:sub
    iV(:,:,j) = eye(len);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Pilot Estimation of the brain response
% First-pass estimation of beta and gamma for each voxel

for v = 1:voxels   

    [beta, gamma, C, AE, sig2e, ksiv, ~, XB, Phi, P] = EstimationVoxel(y(:,:,v), iV, Runs, Init);
  
    % Store results
    Beta(:,v) = beta;
    G(:,v) = gamma;
    Ar(:,v) = AE';
    S2e(v) = sig2e;
    K(:,v) = ksiv;
    CC(:,v) = C;
%    E(:,v) = epsv;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Estimation of the temporal dependence

A = zeros(len,len,parcels);
sig2_eps = zeros(1,parcels);
phi = zeros(p,parcels);

% Estimate parcel-specific matrix Am (AR components)

for m = 1:parcels
    
    ind = (parind == m);
 %   vind = find(ind);
    phi(:,m) = mean(Ar(:,ind));
    sig2_eps(m) = mean(S2e(ind));
     
    Atmp = eye(len);
    for t=1:(len-1)
        Atmp = Atmp + (phi.^t).*(diag(ones(len-t,1),t) + diag(ones(len-t,1),-t));
    end
    A(:,:,m) = (1/(1-phi^2))*Atmp;
    
end


% Estimate rho

rho_l = zeros(nbasis,L,voxels);
sig2_ksi = zeros(L,voxels);
Tksi = zeros(L*nbasis,L*nbasis);

for v=1:voxels

    Ksi2 = zeros(1,L);
    KK = zeros(nbasis,nbasis,L);
    Tksi_l = zeros(nbasis,nbasis,L);
      
    for j=1:N1,

        ksiv = K(:,v);
        ksi = ksiv((j-1)*nbasis*L+1:(j*nbasis*L));
        ksi = reshape(ksi,nbasis,L);

        for l=1:L
            Ksi2(l) = Ksi2(l) + ksi(:,l)'*ksi(:,l);     
            KK(:,:,l) = KK(:,:,l) + ksi(:,l)*ksi(:,l)';
        end
    
    end
    
    for l=1:L

        sig2_ksi(l,v) = Ksi2(l)/(N1*nbasis);
        CVM = KK(:,:,l)/N1;
        
        [~,CM] = cov2corr(CVM);
       
        for k=1:nbasis
            rho_l(k,l,v) = mean(diag(CM,(k-1)));
        end 
        
        Tksi_l(:,:,l) = toeplitz(rho_l(:,l,v));
        Tksi((l-1)*nbasis+1:l*nbasis, (l-1)*nbasis+1:l*nbasis, v) = Tksi_l(:,:,l);
        
    end
    
end


rho = mean(rho_l,3);
Tksi = mean(Tksi,3);
rho_new = rho;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EM-algorithm

% START HERE

maxiter = 10;
tol = 1e-6;
rho_norm = 100;
sig_norm = 100;
niter = 1;
%options = optimset('MaxFunEvals',1000,'Maxiter',1000,'TolX',1e-6,'TolFun',1e-6,'Display','off');
options = optimset('MaxFunEvals',500,'Maxiter',500,'TolX',1e-6,'TolFun',1e-6,'Display','off');

while (rho_norm>=tol && sig_norm>=tol && (niter <= maxiter))

    rho_tmp = zeros(nbasis,L,voxels);
    Tksi_tmp = zeros(L*nbasis,L*nbasis,voxels);

    sig2_ksi_old = sig2_ksi;
    rho_old = rho_new;

     for v=1:voxels

        if (sum(sig2_ksi ==0) == 0) 
            isig2k = 1./sig2_ksi(:,v);
        else
            isig2k = zeros(size(sig2_ksi));
        end
        
        m = parind(v);
        iA = inv(A(:,:,m));
        isig2e = 1/sig2_eps(m);    
        iTksi = pinv(Tksi);
        C = zeros(nbasis,nbasis,L);
        
        k2 = kron(diag(isig2k),eye(nbasis));
        kbg = kron(beta,gamma);
        for j=1:N1

            Bj = inv(XB(:,:,j)'*(isig2e*iA)*XB(:,:,j) + iTksi + k2);
            Rj = eye(len) - Phi(:,:,j)*inv(Phi(:,:,j)'*iV(:,:,j)*Phi(:,:,j))*Phi(:,:,j)'*iV(:,:,j);   
            ksi = Bj*XB(:,:,j)'*(isig2e*iA)*Rj*(y(:,j,v) - XB(:,:,j)*kbg);           
            K((j-1)*L*nbasis+1:j*L*nbasis,v) = ksi;

            for l=1:L
                C(:,:,l) = C(:,:,l) + (ksi((l-1)*nbasis+1:l*nbasis)*ksi((l-1)*nbasis+1:l*nbasis)'); % + Bj((l-1)*nbasis+1:l*nbasis, (l-1)*nbasis+1:l*nbasis)); 
            end

        end

        for l=1:L
            dbstop if error

            rho_tmp(1,l,v) = 1;
            rho_tmp(2:end,l,v) = fminsearch(@minQ,rho_old(2:end,l),options,diag(isig2k),N1,nbasis);

            Tksi_tmp((l-1)*nbasis+1:l*nbasis, (l-1)*nbasis+1:l*nbasis, v) = toeplitz(rho_new(:,l));
            sig2_ksi(l,v) = max(0, trace(pinv(Tksi((l-1)*nbasis+1:l*nbasis, (l-1)*nbasis+1:l*nbasis))*C(:,:,l))/(N1*nbasis));
        end

    end

    rho_new = mean(rho_tmp,3);
    Tksi = mean(Tksi_tmp,3);   

    rho_norm = norm(rho_new - rho_old);
    sig_norm = norm(sig2_ksi - sig2_ksi_old);
    niter = niter + 1;

end

% END HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Estimation of between-subject variance


AA = zeros(L,L);
for l1=1:L
    for l2=1:L
        for j=1:N1
            XBj1 = XB(:,(l1-1)*nbasis+1:l1*nbasis,j);
            XBj2 = XB(:,(l2-1)*nbasis+1:l2*nbasis,j);
            AA(l1,l2) = AA(l1,l2) + trace(XBj1*Tksi((l1-1)*nbasis+1:l1*nbasis, (l1-1)*nbasis+1:l1*nbasis)*XBj1'*XBj2*Tksi((l2-1)*nbasis+1:l2*nbasis, (l2-1)*nbasis+1:l2*nbasis)*XBj2');
        end
    end
end
AA = (AA+AA')/2;

for v =1:voxels

    bb = zeros(L,1);
    m = parind(v);
    
    for j=1:N1
        
        Rj = iV(:,:,j)*(eye(len) - Phi(:,:,j)*inv(Phi(:,:,j)'*iV(:,:,j)*Phi(:,:,j))*Phi(:,:,j)');           
        r_j = Rj*y(:,j) - Rj*XB(:,:,j)*kron(Beta(:,v),G(:,v));
   
        for l=1:L
            XBj = XB(:,(l-1)*nbasis+1:l*nbasis,j);
            bb(l) = bb(l) + (r_j'*XBj*Tksi((l-1)*nbasis+1:l*nbasis, (l-1)*nbasis+1:l*nbasis)*XBj'*r_j - sig2_eps(m)*trace(Tksi((l-1)*nbasis+1:l*nbasis, (l-1)*nbasis+1:l*nbasis)*XBj'*A(:,:,m)*XBj));
        end
    end
    sig2_ksi(:,v) = quadprog(AA,bb,[],[],[],[],zeros(L,1));

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second-pass estimation of beta and gamma for each voxel

% sig2_ksi = zeros(1,voxels);

for v = 1:voxels   

    % Compute inverse of total covariance matrix V
    iV = zeros(len,len,sub);
    D = diag(sig2_ksi(:,v));
    for j=1:sub
        V = XB(:,:,j)*(kron(D,eye(nbasis)))*Tksi*XB(:,:,j)' +  sig2_eps(parind(v))*A(:,:,m); 
        iV(:,:,j) = inv(V);
    end
    
    % Second-pass estimate of beta and gamma
%    [beta, gamma, C, AE, sig2e, ksiv, epsv, XB, Phi] = EstimationVoxel(y(:,:,v), iV, Runs, Init);
    [beta, gamma, C] = EstimationVoxel(y(:,:,v), iV, Runs, Init);
    Beta(:,v) = beta;
    G(:,v) = gamma;
    CC(:,v) = C;
    K(:,v) = ksiv;
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference

sbeta = zeros(L,voxels);
sgamma = zeros(nbasis,voxels);
HT = zeros(1,voxels);

sbeta2 = zeros(L,voxels);
sbeta3 = zeros(L,voxels);


Beta2 = Beta;
whmax = zeros(size(Beta));

for v = 1:voxels   

    M=0;
    Ny =0;
    for j=1:sub
        Rj = eye(len) - Phi(:,:,j)*inv(Phi(:,:,j)'*iV(:,:,j)*Phi(:,:,j))*Phi(:,:,j)'*iV(:,:,j);   
        M = M + XB(:,:,j)'*Rj'*iV(:,:,j)*Rj*XB(:,:,j);
        Ny = Ny + XB(:,:,j)'*iV(:,:,j)*Rj*y(:,j);
    end

    ZB = kron(eye(L),G(:,v));
    VB = ZB'*M*ZB;
    sbeta(:,v) = sqrtm(VB)*Beta(:,v);
    
    for l=1:L
        hh = B*G;
        tmp = hh*Beta(l,v);
        [a,b] = max(abs(tmp));
        Beta2(l,v) = tmp(b);
        whmax(l,v) = b;
    end
    
    sbeta2(:,v) = sqrtm(VB)*Beta2(:,v);

    ZG = kron(Beta(:,v),eye(nbasis));
    VG = ZG'*M*ZG;
    sgamma(:,v) = sqrtm(VG)*(G(:,v)-G0);

    HT(v) = (G(:,v)-G0)'*sqrtm(VG)*(G(:,v)-G0);


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bookkeeping
Res = {};
Res.beta = Beta;        % Voxel-wise beta
Res.gamma = G;          % Voxel-wise gamma
Res.H = B*G;            % Voxel-wise HRF estimates
% Res.Err = Err;
Res.AR = phi;
Res.S2e = sig2_eps;
% Res.K = K;
Res.S2k = sig2_ksi;
Res.rho = rho_new;

Res.sbeta = sbeta;
Res.beta2 = Beta2;

Res.sgamma = sgamma;
Res.sbeta2 = sbeta2;

Res.whmax = whmax;
Res.HT = HT;
Res.K = K;

% Res.sgamma2 = sgamma2;
% Res.SigmaE = SigmaE;
% Res.sPhi = sPhi;
% Res.Vbeta = Vb;
% Res.Vgamma = Vg;


end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SUBFUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Q=minQ(rho,D, N1,nbasis)

rho = [1; rho];

Tksi = toeplitz(rho);

Q = (N1*nbasis*log(det(D)) +N1*log(det(Tksi)) + N1);

end


function [beta, gamma, C, AR, sig2e, ksiv, epsv, XB, Phi, P] = EstimationVoxel(y, iV, Runs, Init)

% Initial values
nbasis = Init.nbasis;
%norder = Init.norder;
Tlen = Init.Tlen;
len = Init.len;
L = Init.L;
N1 = Init.N1;
q = Init.q;
p = Init.p;
lambda = Init.lambda;
hphi = Init.hphi;
B = Init.B;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create model matrices


% Create design matrix

XB = zeros(len, nbasis*L , N1);
for j=1:N1
    Xj = [];
    for l=1:L
        Xjl = tor_make_deconv_mtx3(Runs(:,l,j),Tlen,1);
        Xj = [Xj Xjl(:,1:Tlen)*B];
    end

    XB(:,:,j) = Xj;
end

% Create nuisance matrix
Phi_ij = zeros(len,q);
Phi_ij(:,1) = 1;
Phi_ij(:,2) = (1:len)./len;
% Phi_ij(:,3) = Phi_ij(:,2).^2;
Phi_ij = orth(Phi_ij);      % Orthogonalize matrix

Phi = zeros(len, q, N1);
for j=1:N1
    Phi(:,:,j) = Phi_ij;
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation algorithm


% Step 1: Compute penalty

P = eye(nbasis);

numh = size(hphi,1);
for q = 1:numh,
    P = P - hphi(q,:)'*hphi(q,:);
end


% Step 2: Estimate initial value of gamma

sXRX = 0;
sXRy = 0;

for j=1:N1,         
    Rj = iV(:,:,j)*(eye(len) - Phi(:,:,j)*inv(Phi(:,:,j)'*iV(:,:,j)*Phi(:,:,j))*Phi(:,:,j)'*iV(:,:,j));           
    sXRX = sXRX + XB(:,:,j)'*Rj*XB(:,:,j);
    sXRy = sXRy + XB(:,:,j)'*Rj*y(:,j);
end

h0 = inv(sXRX + N1*lambda*kron(eye(L),P))*sXRy;
h0 = reshape(h0,nbasis,L);
tmp = sum(h0,2);
gamma = tmp./norm(tmp,2);

% Step 3: Iteratively estimate beta and gamma

betaold = zeros(L,1);
beta = ones(L,1);
gammaold = zeros(size(gamma));


cnt = 0;
while ((norm(gamma - gammaold) > 0.0001 || norm(beta-betaold) > 0.0001) && cnt < 1000)
    betaold = beta;
    gammaold = gamma;
    cnt = cnt + 1;

    Z = kron(eye(L),gamma);
    W1 = Z'*sXRX*Z;
    W2 = Z'*sXRy;
    beta = inv(W1'*W1)*W1'*W2;
    
    
    Z = kron(beta,eye(nbasis));
    W1 = Z'*sXRX*Z + N1*lambda*P;
    W2 = Z'*sXRy;

    [U,D,~] = svd(W1);
    lam = diag(D);
   
    [C,~,flag] = fzero(@(C) sum(((W2'*U).^2)'./((lam+C).^2)) - 1, 0);
    if (flag == -6)
        C = 0.1;
    end
        
    W1 = W1 + C*eye(nbasis);
    gamma = inv(W1'*W1)*W1'*W2;
    gamma = gamma./norm(gamma,2);
 
    if (sum(abs(beta)) < 0.001), cnt = 1000; end     % End if all beta are close to zeros
end

% Step 4: Estimation of temporal dependence


AR = zeros(N1,p);
sig2e = zeros(N1,1);
ksiv = zeros(N1*nbasis*L,1);
epsv = zeros(len,1);

for j=1:N1,
        
    Rj = iV(:,:,j)*(eye(len) - Phi(:,:,j)*inv(Phi(:,:,j)'*iV(:,:,j)*Phi(:,:,j))*Phi(:,:,j)');           
    r_j = Rj*y(:,j) - Rj*XB(:,:,j)*kron(beta,gamma);
    ksi = inv(XB(:,:,j)'*XB(:,:,j) + 1000*eye(L*nbasis))*XB(:,:,j)'*r_j;     
    
    epsilon = r_j; 
    epsv = epsv + epsilon;
    
    [a,e] = aryule(epsilon,p);
    
    AR(j,:) = -a(2:end);
    sig2e(j) = e;
    
    ksiv((j-1)*nbasis*L+1:(j*nbasis*L)) = ksi;

end

AR = mean(AR);
sig2e = max(mean(sig2e),0);

end
