function [s,it]=display(a,tab,it)
 
  global display_tree_depth;
  if isempty(display_tree_depth) display_tree_depth=100; end;
  global display_tree_tab_step;
  if isempty(display_tree_tab_step) display_tree_tab_step=3; end;
  global display_tree_show_brackets;
  if isempty(display_tree_show_brackets) display_tree_show_brackets=1; end;

  % ----- tabulation stuff ----------
  print=0;  if nargin==1  tab=0; print=1; it=[]; end;  
  tab_step=display_tree_tab_step;
  % ------ display tree ------------
  [s,it]=display(a.algorithm,tab,it);
  if  display_tree_show_brackets 
    s=[s ones(1,tab+min(2,tab_step))*32 '{\n' ];  
  end
  tab=tab+tab_step;  
  i=length(it); v=it{i}; v=[v 1]; it{i}=v; %% add 1 to last row 
  for i=1:length(a.child)
     [st it]=display(a.child{i},tab,it);
     if length(it{length(it)})<display_tree_depth
      s =[s ones(1,tab)*32 st];       
     end
    
     j=length(it); v=it{j}; v(length(v))=v(length(v))+1; 
     it{j}=v; %% increment last row
  end
  if  display_tree_show_brackets
    if  length(it{length(it)})>=display_tree_depth
     s=[s ones(1,tab)*32 '...\n'];      
    end
    s =[s ones(1,tab-tab_step+min(2,tab_step))*32 '}\n' ];
  end
  tab=tab-tab_step;  
  %% pop index of last item
  i=length(it); v=it{length(it)}; v=v(1:length(v)-1); it{i}=v; 
  if print
   s=s(1:length(s)-2); % remove last carriage return
   it=deal(it(1:length(it)-1)); % remove last index to tree
   disp(sprintf(char(s)));
  end
 
