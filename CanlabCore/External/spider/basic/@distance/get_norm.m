function [val]=get_norm(k,d,ind)
  
% [val]=get_norm(distance,data,[index])
%  Returns the norm of the vectors via the data and
%   distnel specified, optionally one can calculate this only
%   for some of the examples using index
  
if k.calc_on_output d.X=d.Y; end  %% output distnels

if nargin==2
  tmp = get_dim( d);
  ind = [ 1:tmp( 1)];
end
if strcmp( k.dist, 'linear'),
  tmp = get_x( d, ind);  
  val = sqrt( sum( tmp.^2, 2));
else
  %slow implementation but it works
  if length(ind)>3000
   for i=1:length(ind),
  	  val(i) = calc(k.child,get(d,ind(i)));
   end;
  else
   val = calc( k.child, get( d, ind)); val=diag( val)';
   val = sqrt( val');
end
end;



