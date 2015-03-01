function d =  testing( a, d)
%
% 
% Returns the closest codebook vector (or its index) for each given example when the kvq object is in 
% mode 'discriminative','shared' or 'standard'. Otherwise it returns the codebook. 
%

if (strcmp(a.mode,'standard') || strcmp(a.mode,'shared') || strcmp(a.mode,'discriminative'))

  N=get_dim(d);
  N2=get_dim(a.keep);

  D=calc(a.child,d,a.keep);
  [mv,I]=min(D);
else % just return the codebook vectors
  d = a.keep;
end

if(a.return_indices==0)
  d=get(a.keep,I);
else
    d=set_x(d,I);
end

