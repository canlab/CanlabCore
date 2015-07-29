function s=get_name(a)
s=[get_name(a.algorithm)];
eval_name

n=a.func; if length(n)>6 n=[n(1:6) '..']; end;
s=[s  '(' n ')'];

if ~isempty(a.p1) s=[s ' p1=' dis(a.p1)]; end
if ~isempty(a.p2) s=[s ' p2=' dis(a.p2)]; end
if ~isempty(a.p3) s=[s ' p3=' dis(a.p3)]; end
if ~isempty(a.p4) s=[s ' p4=' dis(a.p4)]; end


function s=dis(s)
    
  if isa(s,'double')
    s=num2str(s);
  else
    if isa(s,'algorithm')
      s=get_name(s);
    end
  end

