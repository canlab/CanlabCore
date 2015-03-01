function [str,iter]=display(comp,tab,iter)

  
  global display_tree_array_indexing;
  if isempty(display_tree_array_indexing) 
      display_tree_array_indexing=2; 
  end;
  global display_tree_index_at_front;
  if isempty(display_tree_index_at_front) 
      display_tree_index_at_front=1; 
  end;
      
  if nargin==3
   if isempty(iter)
    iter{1}=[];
   else  
      i=length(iter); 
      iter{i+1}=iter{i};  %% <--- new member of tree
   end
   if display_tree_array_indexing==0 
     if isempty(iter{length(iter)})
      index=[]; 
     else
      index=[num2str(length(iter)-1) ':']; 
     end
     if display_tree_index_at_front
      str=[ index get_name(comp) '\n'];
     else
      str=[ get_name(comp) ' (' index(length(index))  ')\n'];     
     end;
   else
     index=[];
     if ~isempty(iter{length(iter)})
      index= num2str(iter{length(iter)}); 

      a=findstr('  ',index);
      a=setdiff([1:length(index)],a); index=index(a); 
      
      if display_tree_array_indexing==2
    	a=find(index==32); 
        if isempty(a) 
            a=0; 
        end;
	    index=index(max(a)+1:length(index));
      else
        index=['[' index ']']; 
      end
      index=[index ':'];
     end;
     if display_tree_index_at_front
       str=[ index  get_name(comp) '\n'];
     else
       str=[ get_name(comp)  ' ' index(1:length(index)-1)  '\n'];
     end;
   end       
 
  else 
   disp([get_name(comp)]);
  end  
