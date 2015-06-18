function [beta,betastar] = igls_orig(y, x, Vy, type)
% function [betastar] = igls(y, x, Vy)
%
% Variance Component Estimation using IGLS/RIGLS
%
% y = d + cx + epsilon  where  epsilon ~ N(0,sigma*Vy)
%
% d ~ N(0, sigma_d) and c ~ N(0, sigma_c)
%
% Calculate d, c, sigma, sigma_d and sigma_c using Maximum Likelihood
% methods (IGLS) and Restricted Maximum Likelihood methods (RIGLS).
%
% y - matrix T x subjects
% x - matrix T x subjects
% Vy - matrix T x T x subjects
% type = 'i' IGLS
% type = 'r' RIGLS
%
% By Martin Lindquist, April 2007
%
% Example:
% 
% len = 200; sub =20;
% x = zeros(len,sub);
% x(11:20,:) = 2;
% x = x+ normrnd(0,0.1,len,sub);
% c = normrnd(0.5,0.1,sub,1);
% d = normrnd(3,0.2,sub,1);
% y=x;
% for i=1:sub, y(:,i) = d(i) + c(i).*x(:,i) + normrnd(0,0.5,len,1); end;
% Vy = zeros(len,len,sub);
% for i=1:sub, Vy(:,:,i) = eye(len,len).*0.5; end;
% 
% [beta, betastar] = igls(y, x, Vy,'i')
%
c1= clock;
[T, sub] = size(y);             % Length of y vector and Number of subjects

one = zeros(T,1)+1;             % Vector of ones
null = zeros(T,1);              % Vector of zeros

epsilon = 0.001;        
num_iter = 5;

len = sub*T;                                    % Total number of observations
z = reshape(y,len,1);                           % Concatenated data
D = [zeros(len,1)+1 reshape(x,len,1)];          % Design matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 1: Find the OLS solution 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta = pinv(D)*z;                           % Beta values
resid = z - D*beta;                         % Residuals

ystar = [];
if (type == 'i'),           % IGLS
     for i=1:sub,
         tmp = vech(resid(((i-1)*T+1):(i*T))*resid(((i-1)*T+1):(i*T))');                 % Find vech of estimated covariance
         ystar = [ystar; tmp];
     end;     
elseif (type == 'r'),       % RIGLS
    for i=1:sub,
        Dtmp = D((((i-1)*T+1):(i*T)),:);
        rtmp = resid(((i-1)*T+1):(i*T));
        rig = rtmp*rtmp' + Dtmp*inv(Dtmp'*Dtmp)*Dtmp';
        tmp = vech(rig);                              % Find vech of estimated covariance 
        ystar = [ystar; tmp];
    end;
end;
    
clear tmp rtmp rig Dtmp    
G = Create_Design_Eq2(y, x, Vy);                      % Create design matrix for variance estimation
betastar = pinv(G)*ystar;                             % Estimate variance components

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 2: Iterate 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cnt = 0;
betastar_old = betastar+10;

iSigma = zeros(len,len);    

while(cnt < num_iter | sum((betastar-betastar_old).^2)> epsilon),

    num = size(G,1)/sub;
    
    for i=1:sub,
        iSigma(((i-1)*T+1):(i*T),((i-1)*T+1):(i*T)) = ivech(G(((i-1)*num+1):(i*num),:)*betastar);
    end;
    
    beta = inv(D'*iSigma*D)*D'*iSigma*z;          % Beta values
    resid = z - D*beta;                           % Residuals

    ystar = [];

    if (type == 'i'),           % IGLS   
         for i=1:sub,
             tmp = vech(resid(((i-1)*T+1):(i*T))*resid(((i-1)*T+1):(i*T))');                 % Find vech of estimated covariance
             ystar = [ystar; tmp];
         end;
    elseif (type == 'r'),       % RIGLS
         for i=1:sub,
            Dtmp = D((((i-1)*T+1):(i*T)),:);
            rtmp = resid(((i-1)*T+1):(i*T));
            rig = rtmp*rtmp' + Dtmp*inv(Dtmp'*iSigma((((i-1)*T+1):(i*T)),(((i-1)*T+1):(i*T)))*Dtmp)*Dtmp';
            tmp = vech(rig);                                                                 % Find vech of estimated covariance
            ystar = [ystar; tmp];
         end;      
    end;

    clear tmp rtmp rig Dtmp    

    betastar_old = betastar;
    betastar = pinv(G)*ystar;                                         % Estimate variance components

    cnt = cnt+1
end;

c2 = clock;
c2 - c1
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunctions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G] = Create_Design_Eq2(y, x, Vy)
% function [G] = Create_Design(y, x, Vy)
%
% Create Design matrix for estimation of variance componets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T, sub] = size(x);
len = T*(T+1)/2;
ONE = zeros(len,1)+1;

G = [];
H = [];

for i=1:sub,
 
    XX = x(:,i)*x(:,i)';

    Gtmp = zeros(len,3);
    Gtmp(:,1) = ONE;
    Gtmp(:,2) = vech(XX);
    Gtmp(:,3) = vech(Vy(:,:,i));    
    G = [G; Gtmp];
     
end;


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = vech(Mat)
% function V = vech(Mat)
%
% Calculate vech for the matrix Mat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = Mat(logical(tril(ones(size(Mat)))));

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mat = ivech(V)
% function Mat = vech(V)
%
% Calculate the "inverse" of the vech function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len = length(V);
dim = -0.5 + sqrt(0.25 + 2*len);
Mat = zeros(dim,dim);
ind=1;

for i=1:dim
    for j=i:dim
        Mat(j,i)=V(ind);
        ind=ind+1;
    end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
