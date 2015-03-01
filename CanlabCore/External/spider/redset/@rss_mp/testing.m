function d  =  testing(a,d)

%K=calc(a.child.child,d,a.Xsv);
%if length(a.alpha)==0
%  Yest = d.Y*0;
%else
%  Yest = K'*a.alpha;
%  for i=1:length(a.b0) 
%      Yest(:,i)=Yest(:,i)+a.b0(i);     
%  end
%end

%% ---To avoid large kernel matrices, we test in batches---

 sz=get_dim(a.Xsv);   %% <---500x500 point are the maximum for one batch
 if sz==0 sz=1; end;
 batch_size=round((500^2)/sz);
 Yest=[];
 for i=[1:batch_size:get_dim(d)]
   take= [i:min(i+batch_size-1,get_dim(d))];
   if ~isempty(a.alpha) 
     kerMaTemp=get_kernel(a.child.child,get(d,take),a.Xsv);
     Yest=[Yest; (kerMaTemp'*a.alpha)];
   else
     Yest=[Yest ; zeros(length(take),size(d.Y,2))];    
   end
 end

% ..and now add bias to each dimension

  for i=1:length(a.b0) 
    Yest(:,i)=Yest(:,i)+a.b0(i);     
  end


  if a.algorithm.use_signed_output==1  
    Yest=sign(Yest);  
  end  
  d=set_x(d,Yest);    
  d=set_name(d,[get_name(d) ' -> ' get_name(a)]);   
   
   
