function results =  testing(a,d)

%disp(['testing ' get_name(a) '.... '])
  
[m,n] = get_dim(d);
Q=length(a.child);
 dtmp=d; 
 Ytmp=-ones(m,Q);
  for i=1:Q,
    dtmp=set_y(dtmp,Ytmp(:,i));
    dtmp=set_name(dtmp,['Machine ' num2str(i)]);
    Yest(:,i) = get_x(test(a.child{i},dtmp));    
  end;
%% If more than one class in the test set then multi-label

%tmp = sum((Yest==1),2);  THIS IS  SUBTLE! Which one to take?
tmp = sum((Ytmp==1),2);
if ~isempty(find(tmp>1)),
  Yest = sign(Yest);
else
  [r,Ytmp2]=max(Yest,[],2);
  Yest = -ones(size(Yest));
  for i=1:length(Ytmp2),
    Yest(i,Ytmp2(i))=1;
  end;
end;

results=set_x(d,Yest);
if ~isempty(find(tmp>1)),
  results=set_name(results,[get_name(d) ' -> ' get_name(a,1) ' [multi-label]']); 
else 
  results=set_name(results,[get_name(d) ' -> ' get_name(a)]); 
end;





