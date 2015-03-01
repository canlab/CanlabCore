function s=get_name(a,detailed)
s=[get_name(a.algorithm)];
if length(a.child)==1,
    s=[s '(' get_name(a.child{1}) ')'];
else
   s=[s '(' get_name(a.child{1}) '...)']; 
end
