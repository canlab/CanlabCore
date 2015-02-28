function s=get_name(a)  
  
s=[get_name(a.algorithm)];  
if ~strcmp(a.group,'separate')
  s=[s ' [' a.group ']'];
end
