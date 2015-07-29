function K = poly(kern,dat1,dat2,ind1,ind2,kerParam),

% K = poly(d1,d2,ind1,ind2,param), compute the kernel matrix between d1 and d2
%  for a polynomial kernel (<x,z>+1)^param where x is from d1 and z from d2

X2 = get_x(dat2,ind2);
X1 = get_x(dat1,ind1);

% X1 = [X1,X1(:,1)-X1(:,3)];
% X2 = [X2,X2(:,1)-X2(:,3)];


K=(X2*X1'+1).^kerParam;  
