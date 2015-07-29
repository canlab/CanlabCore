% clear all;
t=toy({'l=1000','n=10'})
d=gen(t);
d.Y=exp(d.X(:,1));

% d.X=[1,0,1;1,2,3];
% d.Y=[1;-1];
s=svr(kernel('rbf',2));s.C=Inf;s.epsilon=1e-3;
% s.alpha_cutoff=1e-3;
s.optimizer='libsvm';s.nu=0;

%mex libsvm_regressor_spider.cpp svm.cpp
tic
[r0,sneu]=train(cv(s),d)
toc
loss(r0,'quadratic_loss')
% 
%  s.nu=0;
%  s.optimizer='andre';
%  [r,s0]=train(s,d)
%  loss(r,'quadratic_loss')
% 
