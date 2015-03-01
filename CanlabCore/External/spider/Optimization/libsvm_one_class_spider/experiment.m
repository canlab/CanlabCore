N=500;
c=[5,0;-5,0]
X=1*randn(N,2);
X(1:N/2,:)=X(1:N/2,:)+repmat(c(1,:),N/2,1);
X(N/2+1:end,:)=X(N/2+1:end,:)+repmat(c(2,:),N/2,1);
d=data(X);
% pause
% 
%  [r a]=train(param(one_class_svm({kernel('rbf',1),'optimizer="libsvm"'}),'nu', linspace(0.1,0.9,4)),d);
%  subplot(411); plot(a{1},d);
%  subplot(412); plot(a{3},d);
%  subplot(413); plot(a{3},d);
%  subplot(414); plot(a{4},d);
%  colormap gray;
 
  [r a]=train(one_class_svm({kernel('rbf',1),'optimizer="libsvm"','nu=0.95'}),d);
  [r b]=train(one_class_svm({kernel('rbf',1),'optimizer="libsvm"','nu=0.5'}),d);

  figure(1);plot(a,d);
  figure(2);plot(b,d);