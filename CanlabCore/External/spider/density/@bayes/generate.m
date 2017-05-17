
function d =  generate(a)

 if ~isa(a.child,'cell') r=[]; r{1}=a.child; a.child=r; end;

 ps=a.prior;
 for j=1:length(ps)  
   ps(j)=sum(ps(1:j));
   a.child{j}.l=1;
 end
 ps=[0 ps];
 
 k=length(a.child);
 
 xs=[]; ys=[];
 for i=1:a.l 
   p=rand; j=max(find(p>ps)); %% choose class
   x=get_x(test(a.child{j})); %% generate data with that class
   y=-ones(1,k); y(j)=1;
   xs=[xs ; x]; ys=[ys ; y]; 
 end

 if size(ys,2)==2 ys=ys(:,1); end;

 d=data(get_name(a),xs,ys); 
  
 
