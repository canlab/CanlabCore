function s=get_name(a)
s=[get_name(a.algorithm)];
%if nargin==2 
s=[s '(' get_name(a.child) ')'];
s=[s sprintf(' bags=%d m=%d',a.bags,a.m)];
%end
eval_name
