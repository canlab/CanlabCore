function [str,iter]=display(algo,tab,iter)
 
   
  global display_tree_depth;
  if isempty(display_tree_depth) 
      display_tree_depth=100; 
  end;
  global display_tree_tab_step;
  if isempty(display_tree_tab_step) 
      display_tree_tab_step=3; 
  end;
  global display_tree_show_brackets;
  if isempty(display_tree_show_brackets) 
      display_tree_show_brackets=1; 
  end;



  % <<------- tabulation ---------->>
  prnt=0;  
  if nargin==1  
      tab=0; 
      prnt=1; 
      iter=[]; 
  end;  
  tab_step=display_tree_tab_step;
  
  
  % <<------ tree displaying --------->>
  [str,iter]=display_simple(algo,tab,iter);
  
  if  display_tree_show_brackets 
    str=[str ones(1,tab+min(2,tab_step))*32 '{\n' ];  
  end
  tab=tab+tab_step;  
  
  i1=length(iter); 
  v=iter{i1}; 
  v=[v 1]; %% <--- add 1 to last row 
  iter{i1}=v; %% <--- add 1 to last row 
  for i=1:length(algo.child)
    
    if algo.store_all 
      str=[str ones(1,tab-1)*32 '[FOLD ' num2str(i)  ']\n'];
    end
    
    [st iter]=display(algo.child{i},tab,iter);
     if length(iter{length(iter)})<display_tree_depth
      str =[str ones(1,tab)*32 st];       
     end
    
     i2=length(iter); 
     v=iter{i2}; 
     v(length(v))=v(length(v))+1; 
     iter{i2}=v; %% <--- last row is incremented
  end
  if  display_tree_show_brackets
    if  length(iter{length(iter)})>=display_tree_depth
     str=[str ones(1,tab)*32 '...\n'];      
    end
    str =[str ones(1,tab-tab_step+min(2,tab_step))*32 '}\n' ];
  end
  tab=tab-tab_step;  
  %% <<--- poping the index of last item --->>
  i1=length(iter); 
  v=iter{length(iter)}; 
  v=v(1:length(v)-1); 
  iter{i1}=v; 
  if prnt
   str=str(1:length(str)-2); % <--- remove last line break
   iter=deal(iter(1:length(iter)-1)); % <--- remove last index to tree
   disp(sprintf(char(str)));
  end
 



