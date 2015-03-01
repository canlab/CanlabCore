%demonstration of balanced_ridge
X1=0.3*randn(500,2)+repmat([-1,0],500,1);
X2=randn(20,2)+repmat([1,0],20,1);


d=data([X1;X2],[ones(500,1);-ones(20,1)]);

a0=svm;
a0.C=100;
a0.optimizer='libsvm';

for br=[0.01:1:5];
    clf;
    a0.balanced_ridge=br;
    [r,a]=train(a0,d);
    plot(a);
    pause(0.2);
end