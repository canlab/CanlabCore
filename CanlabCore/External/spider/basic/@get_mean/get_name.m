function s=get_name(a)  
  
s=[get_name(a.algorithm)];  
  
if isa(a.child,'algorithm')
 s=[s '(' get_name(a.child) ')'];  
end  
  
eval_name  
