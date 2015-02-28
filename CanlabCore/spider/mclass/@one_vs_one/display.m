function [str,it]=display(a,tabul,it)


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


  % ----- tabulation stuff ----------
  print=0;  
  if nargin==1  
      tabul=0; 
      print=1; 
      it=[]; 
  end;  
  tab_step=display_tree_tab_step;
 
  % ------ display tree ------------
  [str,it]=display(a.algorithm,tabul,it);
  if  display_tree_show_brackets 
    str=[str ones(1,tabul+min(2,tab_step))*32 '{\n' ];  
  end
  tabul=tabul+tab_step;  
  len=length(it); 
  currIt=it{len}; 
  currIt=[currIt 1]; 
  it{len}=currIt; %% add 1 to last row 
  for i=1:length(a.child)
     [st it]=display(a.child{i},tabul,it);
     if length(it{length(it)})<display_tree_depth
      str =[str ones(1,tabul)*32 st];       
     end
    
     j=length(it); 
     currIt=it{j}; 
     currIt(length(currIt))=currIt(length(currIt))+1; 
     it{j}=currIt; %% increment last row
  end
  if  display_tree_show_brackets
    if  length(it{length(it)})>=display_tree_depth
     str=[str ones(1,tabul)*32 '...\n'];      
    end
    str =[str ones(1,tabul-tab_step+min(2,tab_step))*32 '}\n' ];
  end
  tabul=tabul-tab_step;  
  %% pop index of last item
  len=length(it); 
  currIt=it{length(it)}; 
  currIt=currIt(1:length(currIt)-1); 
  it{len}=currIt; 
  if print
   str=str(1:length(str)-2); % remove last carriage return
   it=deal(it(1:length(it)-1)); % remove last index to tree
   disp(sprintf(char(str)));
  end
 
