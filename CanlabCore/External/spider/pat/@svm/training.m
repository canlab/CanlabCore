function [retDat,algo] =  training(algo,retDat)

if algo.algorithm.verbosity>0
    disp(['training ' get_name(algo) '.... '])
end

opt=algo.optimizer;
if strcmp(algo.optimizer,'default')
    len=length(retDat.Y);
    if len<=200  %% <---if there are more then 200 eamples use svm_light
        opt='andre';
    else
        opt='libsvm';
    end
end
if algo.nu~=0,
    opt='libsvm';
end;


switch opt
    case {'svmtorch'}
      multi = 0;    
      regression = 0;   
      degree = 1;
      gamma = 1; 
      eps = 0.001;
      if strcmp(algo.child.ker,'linear')
       kernelType = 0;
      elseif strcmp(algo.child.ker,'poly')
       kernelType = 1;
       degree = algo.child.kerparam;
      elseif strcmp(algo.child.ker,'rbf')
       kernelType = 2;    
       gamma = algo.child.kerparam;
      end;
      [x y] = get_xy(retDat);   
      C = min(algo.C,10000);
      [alpha, bias0, xSV] = SVMTorch(x,y,regression,multi,kernelType,degree,gamma,C,eps);    
    
    %%<<----------------Andre optimizer-------------------->>
    case {'andre'}   
        
        [KerMa,algo.child]=train(algo.child,retDat); 
        KerMa=KerMa.X;    %% calc kernel
        KerMa=add_ridge(KerMa,algo,retDat); 
        y=get_y(retDat); 
        KerMa=KerMa.*(y*y');
        if( algo.nob == 0)
            [alpha,bias] = quadsolve(KerMa,-ones(size(KerMa,1),1),y',0,algo.C); 
            bias0 = -bias;
        else
            alpha = quadsolve(KerMa,-ones(size(KerMa,1),1),[],0,algo.C); 
        end
        
        alpha= alpha .* y;
        
        
    case {'andre_nob'}   
        
        [KerMa,algo.child]=train(algo.child,retDat); 
        KerMa=KerMa.X;    %% calc kernel
        KerMa=add_ridge(KerMa,algo,retDat); 
        y=get_y(retDat); 
        KerMa=KerMa.*(y*y');
        alpha = quadsolve(KerMa,-ones(size(KerMa,1),1),[],0,algo.C); 
        bias0=0; 
        alpha= alpha .* y;
        %%<<----------------quadprog optimizer-------------------->> 
    case {'quadprog'} 
        
        [KerMa,algo.child]=train(algo.child,retDat); 
        KerMa=KerMa.X;    %% <--- calculate the kernel
        KerMa=add_ridge(KerMa,algo,retDat);  
        len=size(KerMa,1);  
        y=get_y(retDat); 
        KerMa=KerMa.*(y*y');
        opts= optimset('display','off','MaxIter',10000,'LargeScale','off'); 
        [alpha,fVAl,exit,out,lambda] = quadprog(KerMa,- ones(len,1),[],[],y',0, zeros(len,1),algo.C*ones(len,1),[],opts);
        bias0=lambda.eqlin(1);  
        alpha= alpha .* y;
    case {'quadprog_nob'} 
        
        [KerMa,algo.child]=train(algo.child,retDat); 
        KerMa=KerMa.X;    %% <--- calculate the kernel
        KerMa=add_ridge(KerMa,algo,retDat);  
        len=size(KerMa,1);  
        y=get_y(retDat); 
        KerMa=KerMa.*(y*y');
        opts= optimset('display','off','MaxIter',10000,'LargeScale','off'); 
        [alpha,fVAl,exit,out] = quadprog(KerMa,-ones(len,1),[],[],[],[], zeros(len,1),algo.C*ones(len,1),[],opts);
        bias0=0;  
        alpha= alpha .* y;
        %%<<----------------svmlight optimizer-------------------->>  
    case {'svmlight'} 
        
        algo.child.dat=retDat; %% <<-- kernel has to store data now
        [x y]=get_xy(retDat);
        if strcmp(algo.child.ker,'linear') 
            ker=1; 
            param1=1; 
        end;
        if strcmp(algo.child.ker,'poly') 
            ker=2; 
            param1=algo.child.kerparam; 
        end;
        if strcmp(algo.child.ker,'rbf') 
            ker=3; 
            param1=algo.child.kerparam; 
            param1=(2*param1^2);   
        end;
        if strcmp(algo.child.ker,'weighted_linear') 
            ker=1; 
            param1=1; 
            tmp = algo.child.kerparam; 
            x=x .* repmat(tmp,size(x,1),1); %%<--- weight data according to parameters
        end;
        if strcmp(algo.child.ker,'weighted_rbf') 
            ker=4; 
            param1=1; 
            x=get_kernel(algo.child,retDat,retDat); %%<--- calculate kernel
            x = [[1:size(x,1)]' x];
        end;
        if strcmp(algo.child.ker,'custom_fast') | strcmp(algo.child.ker,'custom') | strcmp(algo.child.ker,'from_data')
            x=get_kernel(algo.child,retDat,retDat);  %%<--- calculate kernel
            x = [[1:size(x,1)]' x]; 
            ker=4; 
            param1=1;
        end;
        x=full(x);
        [alphas bias0 ind] = svmlight(x,y,algo.C,algo.ridge,algo.balanced_ridge,ker,param1, ...
            max(0,algo.algorithm.verbosity-1));
        if algo.algorithm.verbosity>1 
            disp('done!'); 
        end;
        alpha=zeros(size(x,1),1); 
        alpha(ind+1)=alphas;

     case {'libsvm'} 

        arglist={};
        x=[];
        y=[];
        svm_type=0;
        kernelType=0;
        degree=3;
        gamma=0;
        coef0=0;
        
        nu=algo.nu;
        if(nu>0)
          svm_type=1;
        end

        s=whos('libsvm_cachesize','global');
        
        if (length(s)>0)
            global libsvm_cachesize;
            cachesize=libsvm_cachesize;
        else
            cachesize=40;
        end
        if algo.algorithm.verbosity>1
         fprintf('Using %d MB Cache for Libsvm\n',cachesize)
        end
        C=algo.C;
        
        
        if(C==Inf)
          C=10^6;
        end
        
        %eps=1e-3;
        %p=0.1;
        %shrinking=1;
        
        eps=1e-5; p=0.001; shrinking=1;
        
        weight_label=[];
        weight=[];
        nr_weight=0;
        
        if strcmp(algo.child.ker,'linear')
          kernelType = 0;
        end;
        if strcmp(algo.child.ker,'poly')
          kernelType = 1; 
          degree = algo.child.kerparam;
          coef0 = 1;
          gamma = 1;
        end;
        if strcmp(algo.child.ker,'rbf'),
          kernelType = 2; 
          sigma = algo.child.kerparam; 
          gamma = 1/(2*sigma^2);
        end;
        
        
        x=get_x(retDat);
        if strcmp(algo.child.ker,'custom'),
          kernelType = 4; 
          K= algo.child.kerparam;
          l = get_dim( retDat);
          x = get_index( retDat);
          x = [ reshape( x, l, 1) [ 1:l]']; % using x to pass indices in Matrix and real indices
        end;
        
        
        y=get_y(retDat); 
        
        arglist={ arglist;};
             
        if( kernelType == 4)
            
          [alpha,xSV,bias0]=libsvm_classifier_spider({'X',x},{'Y',y}, ...
                  {'svm_type',svm_type},{'kerneltype',kernelType},...
                  {'degree',degree},{'gamma',gamma},{'coef0',coef0},...
                  {'cachesize',cachesize},...
                  {'C',C},{'eps',eps},...
                  {'nu',nu},...
                  {'p',p},{'shrinking',shrinking},...
                  {'balanced_ridge',algo.balanced_ridge},{'kmatrix',K});
          algo.Xsv=get(retDat, xSV( :, 2));
        else
          [alpha,xSV,bias0]=libsvm_classifier_spider({'X',x},{'Y',y}, ...
                  {'svm_type',svm_type},{'kerneltype',kernelType},...
                  {'degree',degree},{'gamma',gamma},{'coef0',coef0},...
                  {'cachesize',cachesize},...
                  {'C',C},{'eps',eps},...
                  {'nu',nu},...
                  {'p',p},{'shrinking',shrinking},...
                  {'balanced_ridge',algo.balanced_ridge});
          algo.Xsv=data(xSV);
        end
        
        alpha = alpha * y(1);         
        bias0 = bias0 * y(1); 
        algo.b0 = bias0;
        
        %% code to find which alphas were actually used
        %% in libsvm might be slow but more robust
        %% can switch off using algo.cutoff=-2
        if ~isempty(xSV)
			if algo.alpha_cutoff>-2
	         D=[];
	         for i=1:1000:get_dim(retDat)
	          tak=[i:min(get_dim(retDat),i+999)];
	          D=[D;calc(distance,algo.Xsv,get(retDat,[tak]))];
	         end
             [m1 m2]=min(D);
             if  length(unique(m2))<length(m2)
                 %% SVs are not unique -- there must
                 %%be duplicate data points 
                 for i=1:size(D,2)
                     [m1(i) m2(i)]=min(D(:,i));
                     D(m2(i),:)=10000;
                 end
             end
	         f=m2; 
	         algo.Xsv=retDat;
	         alphas=retDat.Y*0; alphas(f)=alpha; alpha=alphas;
	        end 
		end

        fin = find( abs( alpha)>algo.alpha_cutoff);
        algo.alpha = alpha( fin);
        algo.Xsv = get( algo.Xsv, fin);
        
        if algo.algorithm.do_not_evaluate_training_error
            retDat=set_x(retDat,get_y(retDat)); 
        else
            retDat=test(algo,retDat);
        end
        
        return
end

algo.b0=bias0;
fin=find(abs(alpha)>algo.alpha_cutoff);
algo.alpha=alpha(fin);
algo.Xsv=get(retDat,fin);

if algo.algorithm.do_not_evaluate_training_error
    retDat=set_x(retDat,get_y(retDat)); 
else
    retDat=test(algo,retDat);
end



%% =========================================
%% helper function for xtrain
%% =========================================
function deltab=correct_b(d)


rs=sign(d.X);
max_n_np= length(find( rs==1 &  d.Y==-1));
max_n_pn= length(find( rs==-1 &  d.Y==1));


[rx,I]=sort(-d.X);
rx=-rx;


res=[];

X=d.X;
Y=d.Y;
for i=1:length(rx)
    b=rx(i);
    b0=rx(i);
    r=sign(rx+b0);
    n_np= length(find( r==1 &  Y==-1));
    n_pn= length(find( r==-1 & Y==1));
    res=[res;n_pn,n_np,b0];
end
[a,b]=min( abs(res(:,1)-res(:,2)));
deltab=res(b ,3);

%r=max(res); 
%res=res./repmat(r,length(res),1);
