function r=group2cell(d)
  
%  group2cell  - convert group objects to a flat cell array of groups
%
%  Example: a=group({group({svm svm}) knn}); group2cell(a)
  
  if ~isa(d,'group') r=d; return; end
  
  r=[]; 
  for i=1:length(d.child)
    res=group2cell(d.child{i});
    res=make_cell(res);
    if ~isempty(r) r={r{1:length(r)} res{1:length(res)}};else r=res; end
 end
  
