function s=get_name(a)

s=[get_name(a.algorithm)];
eval_name

if a.seed~=-1       s=[s ' seed=' num2str(a.seed)]; end;
%if a.l~=50          s=[s ' l=' num2str(a.l)];       end;
s=[s ' l=' num2str(a.l)];      
s=[s ' err=' num2str(a.noiselevel)];      
