
%% used to print out parameters  -- needed for most objects! 

 a=struct(a);
 global display_tree_show_defaults; 
 if isempty(display_tree_show_defaults) display_tree_show_defaults=0; end;
 show=display_tree_show_defaults; 

  if isfield(a,'child')&isobject(a.child),
   if isfield(struct(a.child),'ker'),
     s=[s ' ' get_name(a.child)];      
   end
 end;
 if isfield(a,'epsilon')
  if show | a.epsilon~=0
   s=[s ' epsilon=' num2str(a.epsilon)];
  end
 end
 if isfield(a,'C')
  if show | a.C~=Inf
   s=[s ' C=' num2str(a.C)];
  end
 end
 if isfield(a,'ridge')
  if show | a.ridge>1e-13
   s=[s ' ridge=' num2str(a.ridge)];
  end
 end
 if isfield(a,'balanced_ridge')
  if show | a.balanced_ridge~=0
   s=[s ' bal_ridge=' num2str(a.balanced_ridge)];
  end
 end
 if isfield(a,'nu')
  if show | a.nu~=0
   s=[s ' nu=' num2str(a.nu)];
  end
 end
 if isfield(a,'optimizer')
  if show | ~strcmp(a.optimizer,'default')
   s=[s ' optimizer=' a.optimizer];
  end
 end
 if isfield(a,'feat')
  if  ~isempty(a.feat)
      s=[s ' feat=' num2str(a.feat)]; 
  end
 end
 if isfield(a,'speed')
  if show | a.speed
     s=[s ' slow_when<' num2str(a.speed)];
  end
 end 
 if isfield(a,'output_rank')
  if show | a.output_rank==1
     s=[s ' output_rank=' num2str(a.output_rank)];
  end
 end



