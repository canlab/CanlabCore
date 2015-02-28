function s=get_name(a)
s=['distance ']; 
if isempty(a.dist)
  s=[s a.child.ker];
else
  s=[s a.dist];
end

wrt=0;

if iscell(a.distparam),
  s = [s ' ('];
  for i=1:length(a.distparam),
    tmp = a.distparam{i};
    if isa(tmp,'char'),
      s = [s, a.distparam{i} ', '];
      wrt = 1;
    elseif isobject( tmp)
      s = [ s, get_name( tmp) ', '];
      wrt = 1;
    elseif length(tmp)==1&~iscell(tmp),
      s = [s,num2str(tmp) ', '];
      wrt = 1;
    end
  end
  if wrt,   
    s = [s(1:length(s)-2) ') '];
  else
    s = s(1:length(s)-2);
  end
elseif isobject( a.distparam)
  s = [ s ' (' get_name( a.distparam) ')'];
else 
  if ~isempty(a.distparam) 
    if length(a.distparam)==1
      s = [s '=' num2str(a.distparam)  ];
    end
  end
end



