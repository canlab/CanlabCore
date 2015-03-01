function [alpha, b0] = qp_learn(K,Y,C)

  a=qpc(K,Y,-C*(Y==-1),C*(Y==1),zeros(length(K),1),1);
  a=a.*Y;
  alpha=a;

  H = K.*(Y*Y');

  class1 = find (a > max(a)/1e3 & a < 0.99*C & Y == 1);
  class2 = find (a > max(a)/1e3 & a < 0.99*C & Y == -1);
  nAsv = length (class1);
  nBsv = length (class2);
  if (nAsv == 0) & (nBsv==0) 
  b0=(min(H(find(a < 0.99*C & Y== 1),:)*a)-      min(H(find(a < 0.99*C & Y==-1),:)*a))/2;
    warning ('Threshold might be inacurate');
 else

  A = sum(H(class1,:)*a)-nAsv;
  B = sum(H(class2,:)*a)-nBsv;
  b0 = -(A-B)/(nAsv+nBsv);
end

