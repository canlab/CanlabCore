function [res,a] =  training(a,d)
  
 disp(['training ' get_name(a) '.... '])
   
 [m,n,Q] = get_dim(d);

 dtmp=d; 
 Ytmp=get_y(dtmp);
 cache=0;
 if cache
 %% take the kernel of the first child
     koldker =a.child{1}.ker;
     koldkerparam = a.child{1}.kerparam;
     kernelold=kernel(koldker,koldkerparam);
 %% Check if all the children are identical
 for i=1:length(a.child),
     if (isa(a.child{i},'svm')|isa(a.child{i},'stab')|isa(a.child{i},'knn'))
         knewker =a.child{i}.ker;
         knewkerparam = a.child{i}.kerparam;
         kernelnew=kernel(koldker,koldkerparam);
         if ~(kernelold==kernelnew),
             cache=0;
             break;
         end;
     else
         cache=0;
         break;
     end;
 end;
 end
 
 %% if only one type of machine and the kernel
 %% matrix can be stored directly  

 
 if cache==1&(m<1000),
     koldker =a.child{1}.ker;
     koldkerparam = a.child{1}.kerparam;
     kernelold=kernel(koldker,koldkerparam);
     ind = get_index(d);
     Ktmp(ind,ind) = get_kernel(kernelold,d,d);
     a.child{1}.child = kernel('custom',Ktmp);
 end;
 lent = length(a.child);
 
for i=1:Q,
  if (i>lent),
    a.child{i} = a.child{mod(i,lent)+1};
  end;
  dtmp=set_y(dtmp,Ytmp(:,i));
  dtmp=set_name(dtmp,['Machine ' num2str(i)]);
  if i>1&cache,
      a.child{i}.child =kernel('custom',Ktmp);
  end;
  [r{i},a.child{i}]=train(a.child{i},dtmp);
  
  if i>1&cache,
  a.child{i}.child = kernelold;
  end;
  
end;
  
  if cache,
  a.child{1}.child = kernelold;
  end
  
  res=test(a,d);
