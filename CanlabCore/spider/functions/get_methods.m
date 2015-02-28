function [methods,objects,tokens,unique,appear,ind_methods]=get_methods(d)

if ~iscell(d) d=group2cell(d); end;
tokens=[]; objects=[]; methods=[];

appear=sparse(0); unique=[];

for i=1:length(d)

 %% ------------- read in all tokens for given example
 n=d{i}.name;
 f=find(n==' '); f=[0 f length(n)+1];
 t=[];j=1; 
 for k=1:length(f)-1
   tok=n(f(k)+1:f(k+1)-1);  
   if ~strcmp(tok,' ') & ~isempty(tok) %% remove whitespace
     f2=find(tok=='='); if ~isempty(f2) tok=tok(1:min(f2)-1); end;
     t{j}=tok; j=j+1; 
  end;
 end
 
 %% ------------- remove all tokens except object names
 
 k=2; o={t{1}}; yes=0;
 for j=1:length(t)
  if yes  %% this one is an object
    if exist(t{j})==2
      o{k}=t{j}; k=k+1; yes=0;
    end
  end
  if strcmp(t{j},'->')  yes=1; end;  
 end

 %% ------------- find unique tokens and record if they appear
 for j=1:length(o)
   f=strmatch(o{j},unique);
   if isempty(f)  %% ----------------- new unique token found
    l=length(unique); unique{l+1}=o{j};
    appear(i,l+1)=1;
   else %% ------------ old token found.. mark it as appearing
    appear(i,f)=1;
   end
 end

 tokens{i}=t; objects{i}=o; 
end


%% ------- assume real methods are ones which don't occur all the time
occur=sum(appear); 
ind_methods=find(occur~=length(appear));
methods=unique(ind_methods);
if isempty(methods)  %%         OOPS! no methods left: pick only one method
  ind_methods=length(unique);
  methods={unique{ind_methods}}; 
end; 
