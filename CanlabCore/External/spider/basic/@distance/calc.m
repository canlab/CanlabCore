function [KM] = calc(k,d1,d2,ind1,ind2)
 
% [KM] = calc(k,d1,[d2],[ind1],[ind2])
%
% Calculate the distance with specific indexed data from 
% datasets d1 and (optionally) d2.
% If d2 is not given it is assumed to make the matrix between d1 and itself
  
  global caching_kernel;

  if nargin<3, d2=[];end;
  if isempty(d2) d2=d1; end;
  if nargin<4   
    ind1=[1:length(get_index(d1))];
  end
  if nargin<5
    ind2=[1:length(get_index(d2))];
  end

  if k.calc_on_output d1.X=d1.Y; d2.X=d2.Y; end %% output distance
  
  index1 = get_index(d1);
  index2 = get_index(d2);
  if strcmp(k.dist,'custom')
   % the distance matrix is stored in param{1} of distance   
   KM = k.distparam;
   KM = KM(index2(ind2),index1(ind1));
   return;
  end
 if strcmp(k.dist,'custom_fast')
   global K; 
   KM = K(index2(ind2),index1(ind1));
   return;
 end
 if k.kercaching==1 & (~isa(d1,'data_global')|~isa(d2,'data_global')),
   error('Caching distance only for data_global');
 else
   if (k.kercaching==1)&~ismyX(d1)&~ismyX(d2), 
     if caching_kernel==k,
       global K;
       index1 = get_index(d1);
       index2 = get_index(d2);
       KM = K(index2(ind2),index1(ind1)); 
       return;
     else
       global K;
       %global X;global Y;
       dtmp=data_global('get_distance tmp',[],[]);
       [l n q]  = get_dim(dtmp);
       ind1=[1:l]; ind2=[1:l];
       caching_kernel=k;
     end
   end
 end
 
 paramtmp=[];
 if iscell(k.distparam),
   for i=1:length(k.distparam),
     paramtmp{i} = k.distparam{i};
   end
 else
   paramtmp = k.distparam;
 end

 if k.kercaching~=1|ismyX(d1)|ismyX(d2),
   if ~isa(k.child,'kernel'),
     KM = feval(k.dist,k,d1,d2,ind1,ind2,paramtmp); 
     KM=real(KM);  
   else
     KM = sqrt((get_norm(k.child,d1,ind1).^2*ones(1,length(ind2)))'+get_norm(k.child,d2,ind2).^2*ones(1,length(ind1))-2*get_kernel(k.child,d1,d2,ind1,ind2));   
     KM=real(KM);
   end;
 else
   if ~isa(k.child,'kernel'),
     K = feval(k.dist,k,dtmp,dtmp,ind1,ind2,paramtmp); 
     K=real(K);  
   else
     K = sqrt(get_norm(k.child,dtmp,ind1).^2*ones(1,length(ind2))'+get_norm(k.child,dtmp,ind2).^2*ones(1,length(ind1))-2*get_kernel(k.child,dtmp,dtmp,ind1,ind2)); 
     K=real(K);
   end  
   KM = K(get_index(d2),get_index(d1));
 end;
 
 if k.output_distance  %% ------------ convert distance to a distance
   error('output distance: not yet implemented')
 end
 