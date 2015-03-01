
function [Rnew,Ki]=finvupdate(a,R,k_newvsold,knew)
k=k_newvsold;
k22=knew;
e=R'*R*k;
g= 1/   sqrt(k22-k'*R'*R*k); 
Rnew=[ R,zeros(length(k),1);-g*e',g];
if(nargout==2)
    Ki=Rnew'*Rnew;
end