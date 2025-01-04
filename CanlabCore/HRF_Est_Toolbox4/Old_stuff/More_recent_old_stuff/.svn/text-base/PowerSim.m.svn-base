N = 1000;
M = 10;
P = zeros(N,M);
B = zeros(N,M);

Yind = zeros(60,1);
Yind(21:40) = 1;
X2 = Yind;
sigma = 0.5;
tstar = tinv(1-0.05,59);

for i=1:M,
    for rep = 1:N,

        X = zeros(60,1);
        X((20+i):(40+i-1)) = 1;
        Y = Yind + normrnd(0,sigma,60,1);
        b = pinv(X)*Y;
        e = Y-X*b;
        sig = sqrt(e'*e/59);
        t = b./(sig/(X'*X));
        delta = 1/sig;
        Pow = 1- nctcdf(tstar,59,delta);

        b2 = pinv(X2)*Y;
        e2 = Y-X2*b2;
        sig2 = sqrt(e2'*e2/59);
        t2 = b2./(sig2/(X2'*X2));
        delta2 = 1/sig2;
        Pow2 = 1- nctcdf(tstar,59,delta2);

       
        P(rep,i) = Pow2-Pow;
        B(rep,i) = b-1;

    end;
end;
