function d =  testing(a,d)

  eval([a.func ';']);
  d=set_name(d,[get_name(d) ' -> ' get_name(a)]);
   
