function D = norm(a,d1,d2,ind1,ind2,kerparam),
   
	if isempty(kerparam)
		kerparam=2;
	end
  	x=get_x(d2,ind2)'; y=get_x(d1,ind1)';   
	[d,m] = size(x);
	[d,n] = size(y);
	Im = ones(1,m);
	In = ones(1,n);
	y = reshape(y,[d 1 n]);
	if kerparam == 1
		D = reshape(sum(abs(y(:,Im,:)-x(:,:,In)),1),[m n]);
	else				
		D = reshape(sum((y(:,Im,:)-x(:,:,In)).^kerparam,1),[m n]);
		D=D.^(1/kerparam);	
	end

% a bit slow but it works..for now
%%   l1=size(x,1); l2=size(y,1); 
%%   D=zeros(l1,l2);
%%    
%%   for i=1:l1
%%     for j=1:l2
%%       D(i,j)=norm(x(i,:)-y(j,:),kerparam);
%%     end
%%   end  
