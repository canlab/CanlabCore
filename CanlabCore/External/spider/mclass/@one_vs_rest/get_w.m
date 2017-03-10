function w = get_w(sv,method) 

if nargin==1 method=1; end;
  
w=[];
for i=1:length(sv.child)  
  wnew= get_w(sv.child{i});  
  wnew=wnew/norm(wnew);
  w=[w;wnew];
end;

w=sum(w);

%  w=max(w);  -- maybe this could be good?



