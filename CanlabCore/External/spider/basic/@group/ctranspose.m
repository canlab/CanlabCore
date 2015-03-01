function d = ctranspose(d)
  
if length(d.child{1}.child)>1
  d2=[]; l=length(d.child);
  for i=1:l
    dd=d.child{i}.child;
    for j=1:length(d.child{i}.child)
      d2{j}{i}=dd{j};
    end    
  end
  d=[];
  for i=1:length(d2)
    d{i}=group(d2{i});
  end
  d=group(d);
end
