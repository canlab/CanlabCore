function r=group2vec(d)
   
%  group2vec   - convert group objects Y output to a flat cell array of outputs
%
%  Example: a=loss(train(cv(svm),gen(toy))); group2vec(a)
  
  
  r=[];
  r2=group2cell(d);
  for i=1:length(r2)
    r=[r ; r2{i}.Y];
  end
