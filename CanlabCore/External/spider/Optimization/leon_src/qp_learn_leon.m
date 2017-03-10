function [alpha, b0] = qp_learn_leon(K,Y,C)
  a=svqp(K,Y,-C*(Y==-1),C*(Y==1),zeros(length(K),1),1);
  a=a.*Y;
  alpha=a;

  H = K.*(Y*Y');

  classAsvi = find (a > max(a)/1e3 & a < 0.99*C & Y == 1);
  classBsvi = find (a > max(a)/1e3 & a < 0.99*C & Y == -1);
  nAsv = length (classAsvi);
  nBsv = length (classBsvi);
  if (nAsv == 0) & (nBsv==0) 
  b0=(min(H(find(a < 0.99*C & Y== 1),:)*a)-      min(H(find(a < 0.99*C & Y==-1),:)*a))/2;
    warning ('Threshold might be inacurate');
 else

  A = sum(H(classAsvi,:)*a)-nAsv;
  B = sum(H(classBsvi,:)*a)-nBsv;
  b0 = -(A-B)/(nAsv+nBsv);
end

