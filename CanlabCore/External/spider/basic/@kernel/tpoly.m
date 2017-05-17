function K = tpoly(k,d1,d2,ind1,ind2,kerparam),

% K = tpoly(d1,d2,ind1,ind2,param), compute the kernel matrix between d1 and d2
%  for a t-polynomial kernel where x is from d1 and z from d2
% Example of definition of tpoly kernel:
% kerp = {[1 2 3] [4 5] [1 6]};
% k = kernel('tpoly',kerp);
%
% then when get_kernel is called, it computes: k(x,y)=<x1,y1><x2,y2>...<xn,yn>
% where n is equal to 3 here and is the number of vectors in kerp (cell)
% and xi (resp. yi) is the vector x restricted to the indices stored in
% the vectors kerp{i}. This kernel is used in nfe.

K=ones(size(get_x(d2,ind2)*get_x(d1,ind1)'));
for i=1:length(kerparam),
    findex = kerparam{i};
    K=(get_x(d2,ind2,findex)*get_x(d1,ind1,findex)'+1).*K;  
end;
