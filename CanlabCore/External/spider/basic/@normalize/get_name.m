function s=get_name(a)

s=[get_name(a.algorithm)];  %% print name of algorithm
s=[s ' scale_type=' num2str(a.scale_type) ];  
if a.sigmoid
 s=[s ' sig=' num2str(a.sig) ];  
end
eval_name                   %% print values of hyperparameters
