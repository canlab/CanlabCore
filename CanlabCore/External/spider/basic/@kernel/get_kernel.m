function [KM] = get_kernel(kern,dat1,dat2,ind1,ind2)

% [KM] = get_kernel(k,d1,[d2],[ind1],[ind2])
%
% Calculate the kernel with specific indexed data from 
% datasets d1 and (optionally) d2.
% If d2 is not given it is assumed to make the matrix between d1 and itself
  
  global caching_kernel;

  if nargin<3, 
      dat2=[];
  end;
  if isempty(dat2) 
      dat2=dat1; 
  end;

  if nargin<4   
    ind1=[1:length(get_index(dat1))];
  end
  if nargin<5
    ind2=[1:length(get_index(dat2))];
  end

  if kern.calc_on_output %% output kernels 
      dat1.X=dat1.Y; 
      dat2.X=dat2.Y; 
  end 
  
  index1 = get_index(dat1);
  index2 = get_index(dat2);
  if strcmp(kern.ker,'custom')
   % the kernel matrix is stored in param{1} of kernel   
   KM = kern.kerparam;
   KM = KM(index2(ind2),index1(ind1));
   return;
  end
 if strcmp(kern.ker,'custom_fast')
   global K; 
   KM = K(index2(ind2),index1(ind1));
   return;
 end
 if kern.kercaching==1 & (~isa(dat1,'data_global')|~isa(dat2,'data_global')),
   error('Caching kernel only for data_global');
 else
   if (kern.kercaching==1)&~ismyX(dat1)&~ismyX(dat2), 
     if caching_kernel==kern,
       global K;
       index1 = get_index(dat1);
       index2 = get_index(dat2);
       KM = K(index2(ind2),index1(ind1)); 
       return;
     else
       global K;
       %global X;global Y;
       dtmp=data_global('get_kernel tmp',[],[]);
       [l n q]  = get_dim(dtmp);
       ind1=[1:l]; ind2=[1:l];
       caching_kernel=kern;
     end
   end
 end
 
  paramTemp=[];
  if iscell(kern.kerparam),
    for i=1:length(kern.kerparam),
      paramTemp{i} = kern.kerparam{i};
    end
  else
    paramTemp = kern.kerparam;
  end
  if kern.kercaching~=1|ismyX(dat1)|ismyX(dat2),
    

    %% --- normal Kernel calculation ----------------------
    KM = feval(kern.ker,kern,dat1,dat2,ind1,ind2,paramTemp);

  else
    K = feval(kern.ker,kern,dtmp,dtmp,ind1,ind2,paramTemp);
    KM = K(get_index(dat2),get_index(dat1));
  end;
  
  if kern.output_distance  %% ------------ convert kernel to a distance
    kern.output_distance=0;
    Kdn = get_norm(kern,dat1,ind1).^2; 
    Kn = get_norm(kern,dat2,ind2).^2;  
    KM = ones(length(Kn),1)*Kdn' + Kn*ones(1,length(Kdn)) - 2*KM;
  end
