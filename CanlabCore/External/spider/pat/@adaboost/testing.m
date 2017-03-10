function d =  testing(a,d)
  
 n= get_dim(d); 
 r = zeros(n,1);
 for i=1:a.kmax 	
 	r =  r + a.alpha(i).*get_x(test(a.child{i},d));
 end

    
 r=sign(r);
%  if(c > 1)
%  	L=-ones(c)+2*eye(c);
% 	for i=1:n
% 		[mi ind]=min(sum(abs(repmat(r(i,:),c,1)-L),2));
% 		r(i,:)=L(ind,:);
% 	end	
%  end
 
 d=set_x(d,r); 
 d=set_name(d,[get_name(d) ' -> ' get_name(a)]); 
