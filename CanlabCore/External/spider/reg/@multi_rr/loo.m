
function [r]=loo(a,d)

% returns the leave one out error in the form of the left-out
% predictions for every point in the training set

res=d;

K=d.X; Y=d.Y; y=[]; 
A=inv(K'*K);
err=[];
for i=1:length(Y)
 t1=Y(i)- K(i,:)*A*K'*Y;
 t2=1- K(i,:)*A*K(i,:)';
 y=[y ; -(t1/t2-Y(i))]; 
% err=[err (t1/t2)^2];
end 

r=set_x(res,y);


