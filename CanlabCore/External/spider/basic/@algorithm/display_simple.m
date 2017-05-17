function [str,iter]=display_simple(comp,tab,iter)
  
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
      a=find(index==32); 
      a=a(1:2:length(a)); 
      a=setdiff([1:length(index)],a); 
      index=index(a); 
      
      if display_tree_array_indexing==2
         index= num2str(iter{length(iter)});
	     fTmp=find(index==' '); 
         fTmp=max(setdiff(fTmp,length(index)));
	        if ~isempty(fTmp)
	            index=index(fTmp+1:length(index));
        	end
      end
      index=['[' index ']:'];
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
