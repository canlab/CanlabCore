%demonstration of custom kernel for libsvm
t=toy({'l=3000','n=20'})
d=gen(t);

%d.X=[1:20]';

k = kernel( 'rbf', 2);
K = calc( k, d);

ind = 3*randperm(10);
% d2= get(d, ind);
 d2= d;
clear d

kc = kernel( 'custom',K);
s = svm(kc);
s.optimizer='libsvm';
[r,a1]=train(s,d2);

s = svm(k);
s.optimizer='libsvm';
[r,a2]=train(s,d2);



[a1.alpha,a2.alpha]
sum( a1.Xsv.X - a2.Xsv.X)
