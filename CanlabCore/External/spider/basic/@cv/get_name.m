function s=get_name(a)
 
global display_tree_show_defaults; 
if isempty(display_tree_show_defaults) 
    display_tree_show_defaults=0; 
end;

show=display_tree_show_defaults; 
s=[get_name(a.algorithm)];
s=[s ' folds=' num2str(a.folds)];
if a.repeats>1
 s=[s ' repeats=' num2str(a.repeats)];
end
if show | a.balanced
 s=[s ' balanced=' num2str(a.balanced) ];
end
eval_name
