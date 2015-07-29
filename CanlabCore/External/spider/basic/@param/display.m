function [s,it]=display(a,tab,it)


  global display_tree_depth;
  if isempty(display_tree_depth) display_tree_depth=100; end;
  global display_tree_tab_step;
  if isempty(display_tree_tab_step) display_tree_tab_step=3; end;
  global display_tree_show_brackets;
  if isempty(display_tree_show_brackets) display_tree_show_brackets=1; end;
 


% ----- tabulation stuff ----------
tab_step=display_tree_tab_step;
print=0;  if nargin==1  tab=0; print=1; it=[]; end;  
  
if ~a.store_all %% simple output
   [s,it]=display_simple(a,tab,it);
else  
  
  [s,it]=display(a.algorithm,tab,it);
  if a.force_no_train
   s=s(1:length(s)-2); % remove last carriage return   
   s=[s ' force_no_train=' num2str(a.force_no_train)  '\n'];
  end
  
  if  display_tree_show_brackets 
    s=[s ones(1,tab+min(2,tab_step))*32 '{\n' ];  
  end
  tab=tab+tab_step;  
  i=length(it); v=it{i}; v=[v 1]; it{i}=v; %% add 1 to last row 
  %% --------  main loop to determine params   -------------
  if isa(a.child,'cell') m=a.child{1}; else m=a.child; end; 
  if ~isa(a.values,'cell') val{1}=a.values; else val=a.values; end
  if ~isa(a.param,'cell') p{1}=a.param; else p=a.param; end 
  for i=1:length(val) sz(i)=length(val{i}); end;
  tot=prod(sz); %% find permutations of hyperparameters  
  for i=1:tot %% all permutations
    vars=num2choice(a,i,sz); %% hyperparams for this iteration      
    for j=1:length(p)      %% set hypers       
      v=(val{j}(vars(j)));
      eval([ 'm.' p{j} '=v;']);   
    end 
      
      if isa(a.child,'cell')
       [st it]=display(a.child{i},tab,it);
       else
       [st it]=display(m,tab,it);
     end
     if length(it{length(it)})<display_tree_depth
      s =[s ones(1,tab)*32 st];       
     end
    
     j=length(it); v=it{j}; v(length(v))=v(length(v))+1; 
     it{j}=v; %% increment last row
  end
  %% ---------------------
  if  display_tree_show_brackets
    if  length(it{length(it)})>=display_tree_depth
     s=[s ones(1,tab)*32 '...\n'];      
    end
    s =[s ones(1,tab-tab_step+min(2,tab_step))*32 '}\n' ];
  end
  tab=tab-tab_step;  
  %% pop index of last item
  i=length(it); v=it{length(it)}; v=v(1:length(v)-1); it{i}=v; 
end
if print
  s=s(1:length(s)-2); % remove last carriage return
  it=deal(it(1:length(it)-1)); % remove last index to tree
  disp(sprintf(char(s)));
end
