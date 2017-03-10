function s=get_name(a)
s=[get_name(a.algorithm)];
%if nargin==2 
s=[s '(' get_name(a.child) ')'];
%end
eval_name
