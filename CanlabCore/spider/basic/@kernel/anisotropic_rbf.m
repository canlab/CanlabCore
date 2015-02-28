function K = anisotropic_rbf(kern,dat1,dat2,ind1,ind2,kerParam),

% K = rbf(d1,d2,ind1,ind2,param), compute the kernel 
%     matrix between d1 and d2
% for a rbf kernel exp(-||x-z||^2/(2*param^2)) 
%     where x is from d1 and z from d2

w = [];
for i = 1:length(kerParam)-1
    w = [w,kerParam{i}];
end
v = kerParam{end};

X2 = get_x(dat2,ind2);
X1 = get_x(dat1,ind1);

for i  =1:size(X1,1)
    for j = 1:size(X2,1)
         s = sum(w.*((X1(i,:) - X2(j,:)).^2));
         K(j,i) = v*exp(-0.5*s);   
     end
end
