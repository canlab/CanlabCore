function [retK] =  add_ridge(retK,algo,dat)
  
  %% <<-----add ridge----->>
  if algo.ridge>0
    retK=retK+sparse(eye(size(retK))*algo.ridge);
  end
    
  %% <<-----add balanced ridge------>>
  if algo.balanced_ridge>0
  % <<-----divide by no. of examples in each----->>
  % Note:
  % It would be possible to calculate the median of diagonal 1,2
  % (dicarded because of implementation difficulties in svm_light)
 
    y=get_y(dat); 
    len=length(y);
    
    rid=diag(retK);
    
%    COULD TAKE MEDIANS
%    fin=find(y==1); rid(fin)=rid(fin)*0+median(rid(fin))*(length(fin)/len)*algo.balanced_ridge;
%    fin=find(y==-1); rid(fin)=rid(fin)*0+median(rid(fin))*(length(fin)/len)*algo.balanced_ridge;
%    OR DO BELOW INSTEAD
    fin=find(y==1);   
    rid(fin)=rid(fin)*0+(length(fin)/len)*algo.balanced_ridge;
    
    fin=find(y==-1); 
    rid(fin)=rid(fin)*0+(length(fin)/len)*algo.balanced_ridge;
    retK=retK+diag(rid);
         
  end
 
