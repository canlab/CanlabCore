function [d,a] =  training(a,d)

disp(['training ' get_name(a) '.... '])

switch a.optimizer

    case {'andre_c'}

        K=train(a.child,d); K=K.X;  % calc kernel
        m=length(K);
        alpha = quadsolve(K,-ones(m,1),[],0,a.C);
        a.alpha=alpha;
        a.b0=0;
        a.Xsv=d;

    case {'quadprog_c'}

        K=train(a.child,d); K=K.X;  % calc kernel
        m=length(K);
        opts=optimset('LargeScale','off');
        [alpha,fval,exitflag,output,lambda]=quadprog(K,-ones(m,1),[],[],[],[],zeros(m,1),ones(m,1)*a.C,[],opts);
        a.alpha=alpha;
        a.b0=0;
        a.Xsv=d;

    case {'quadprog'}

        K=train(a.child,d); K=K.X;  % calc kernel
        m=length(K);
        opts=optimset('LargeScale','off');
        [alpha,fval,exitflag,output,lambda]=quadprog(K,zeros(m,1),[],[],ones(1,m),1,zeros(m,1),ones(m,1)/(a.nu*m),[],opts);
        a.alpha=alpha;
        a.b0=lambda.eqlin;
        a.Xsv=d;

    case {'andre'}
        error('Not implemented yet');


    case {'libsvm'}

        %% default values for libsvm
        cachesize = 40;
        eps=1e-5; p=0.00001; shrinking=1;
        gamma=1;
        C=1;
        svm_type = 2;
        degree = 0;
        coef0 = 0;
        nu = a.nu;
        kertmp = a.child;
        if strcmp(kertmp.ker,'linear')
            kernelType = 0;
        end;
        if strcmp(kertmp.ker,'poly')
            kernelType = 1;
            degree = kertmp.kerparam;
            %deg = optcell.deg;
            coef0 = 1;
        end;
        if strcmp(kertmp.ker,'rbf'),
            kernelType = 2;
            sigma = kertmp.kerparam;
            %sigma = optcell.sigma;
            gamma = 1/(2*sigma^2);
        end;
        if strcmp(kertmp.ker,'custom'),
            kernelType = 4;
            K= a.child.kerparam;
            l = get_dim( retDat);
            x = get_index( retDat);
            x = [ reshape( x, l, 1) [ 1:l]']; % using x to pass indices in Matrix and real indices
        end;

        s=whos('libsvm_cachesize','global');

        if (length(s)>0)
            global libsvm_cachesize;
            cachesize=libsvm_cachesize;
        else
            cachesize=40;
        end
        if a.algorithm.verbosity>1
            fprintf('Using %d MB Cache for Libsvm\n',cachesize)
        end


        Y=get_y(d); 
        X=get_x(d);

        if( kernelType == 4)

            [alpha,xSV,bias0]=libsvm_one_class_spider({'X',X},{'Y',Y}, ...
                {'svm_type',svm_type},{'kerneltype',kernelType},...
                {'degree',degree},{'gamma',gamma},{'coef0',coef0},...
                {'cachesize',cachesize},...
                {'C',C},{'eps',eps},...
                {'nu',nu},...
                {'p',p},{'shrinking',shrinking},...
                {'balanced_ridge',0},{'kmatrix',K});
            a.Xsv=get(retDat, xSV( :, 2));
        else
            [alpha,xSV,bias0]=libsvm_one_class_spider({'X',X},{'Y',Y}, ...
                {'svm_type',svm_type},{'kerneltype',kernelType},...
                {'degree',degree},{'gamma',gamma},{'coef0',coef0},...
                {'cachesize',cachesize},...
                {'C',C},{'eps',eps},...
                {'nu',nu},...
                {'p',p},{'shrinking',shrinking},...
                {'balanced_ridge',0});
            a.Xsv=data(xSV);
        end
        
        b0=-bias0;% from libsvm

        %% code to find which alphas were actually used
        %% in libsvm might be slow but more robust
        %% can switch off using a.cutoff=-2
        if ~isempty(xSV)
            if a.alpha_cutoff>-2
                D=[];
                for i=1:1000:get_dim(d)
                    tak=[i:min(get_dim(d),i+999)];
                    D=[D;calc(distance,a.Xsv,get(d,[tak]))];
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
                a.Xsv=d;
                alphas=d.Y*0; alphas(f)=alpha; alpha=reshape(alphas,length(alphas),1);
            end
        end

        fin = find( abs( alpha)>a.alpha_cutoff);
        a.alpha = alpha( fin);
        a.Xsv = get( a.Xsv, fin);

end
if a.algorithm.do_not_evaluate_training_error==1
    d=set_x(d,get_y(d));
else
    d=test(a,d);
end
