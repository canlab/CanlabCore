TR = 0.5;
len = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute hrf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=0.1:TR:30;
a1 =6;
a2 = 12;
b1 = 0.9;
b2 =0.9;
c = 0.35;

d1 = a1*b1;
d2 = a2*b2;

h = ((t./d1).^a1).*exp(-(t-d1)./b1) - c*((t./d2).^a2).*exp(-(t-d2)./b2);
h = h./max(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X = zeros(len,4);
% X(:,1) = 1;
% X(:,2) = (1:len)/len;
% X(:,3) = X(:,2).^2;
% 
% Run = zeros(40,5);
% Run(1:20,:) = 1;
% Run = reshape(Run,200,1);
% q = conv(Run,h);
% X(:,4) = q(1:len);


X = zeros(len,2);
X(:,1) = 1;
Run = zeros(200,1);
Run(50) = 1;
Run(150) = 1;
q = conv(Run,h);
X(:,2) = q(1:len);


Run = zeros(200,1);
Run(50) = 1;
Run(160) = 1;
tc = conv(Run,h);
tc = tc(1:len) + normrnd(0,0.3,200,1);


beta = pinv(X)*tc;
e = tc-X*beta;
sigma = sqrt((1/(len-size(X,2)))*sum(e.^2));
%c = [0 0 0 1];
c = [0 1];
se = sigma.*sqrt(c*inv(X'*X)*c');
t = beta/se;
pval = 2*tcdf(-abs(t),len-4);
df = len - 4;

[z sres sres_ns] = ResidScan(e, 4);
[b bias pl pc pe] = BiasPowerloss(tc, X,c,beta,df,z,0.05);

beta
b
bias
pl
