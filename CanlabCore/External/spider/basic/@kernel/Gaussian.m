function K = gaussian(kern,d1,d2,ind1,ind2,kerParam),

% K = gaussian(d1,d2,ind1,ind2,param), compute the kernel 
%     matrix between d1 and d2
% for a gaussian kernel exp(-||x-z||^2/(2*param^2)) 
%     where x is from d1 and z from d2

   K=get_x(d2,ind2)*get_x(d1,ind1)';  
   kertmp=kernel;
   Kdn = get_norm(kertmp,d1,ind1).^2; 
   Kn = get_norm(kertmp,d2,ind2).^2;  
   K = ones(length(Kn),1)*Kdn' + Kn*ones(1,length(Kdn)) - 2*K;
   
   [numEx vDim oDim]=get_dim(d1); 
   sigma=kerParam;
   
%   K = 1/(2*(2*pi)^(vDim/2)*sigma) * exp(-0.5*sigma^2*K);

K= 1/(sigma*(2*pi)^(vDim/2))*exp(-K/(2*sigma^2)); 