function s=get_name(a)
s=[get_name(a.algorithm)];

if ~iscell(a.child)
 s=[s '(' get_name(a.child) ')'];
end

if 0
  if ~isa(a.values,'cell') val{1}=a.values; else val=a.values; end
  if ~isa(a.param,'cell') p{1}=a.param; else p=a.param; end 
  for i=1:length(p)
    %index=num2str(val{i}); 
    %t=find(index==32); t=t(1:2:length(t)); 
    %t=setdiff([1:length(index)],t); index=index(t); 
    %s=[s ' ' p{i}  '=[' index ']']; 
    s=[s ' ' p{i}  ]; 
  end
end
  if a.force_no_train
    s=[s ' force_no_train=' num2str(a.force_no_train)];
  end
