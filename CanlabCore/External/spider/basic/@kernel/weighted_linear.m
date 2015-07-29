function K = weighted_linear(k,d1,d2,ind1,ind2,kerparam), 
   
K = get_x(d2,ind2);   K = K .* repmat(kerparam{1},size(K,1),1); 
K2 =get_x(d1,ind1);   K2 = K2 .* repmat(kerparam{1},size(K2,1),1); 
K= K*K2'; 
