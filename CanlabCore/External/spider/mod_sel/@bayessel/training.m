function [d,a] =  training(a,d)


s=inline('1./(1+exp(-100*u))','u');
s1=inline('(1-s(u)).*s(u)','u','s');
s2=inline('(1-s(u)).*s1(u,s) - s1(u,s).*s(u)','u','s','s1');
r=inline('u.*s2(u,s,s1) + 2*s1(u,s)','u','s','s1','s2');

numEx=get_dim(d);
finished=0;
eps=1e-3;
tau=0.1;
s0=a.child;
evidence_max=-inf;
lambda=1;
alpha=1;
beta=1;
lambda_old=lambda;
alpha_old=alpha;
beta_old=beta;
max_iter=20;
i=1;
if isa(s0,'svm')
	while ~finished	
		lambda
		if strcmp(a.type,'L1')
			if a.use_balanced_C ~= 0
				s0.balanced_ridge=lambda;
				s0.optimizer='libsvm';
			else
				s0.C=1/lambda;								
				if isinf(s0.C) | (s0.C > 1e5)
					s0.C=1e5;
				end;
			end;			
			[t s0]=train(s0,d);
			s0.algorithm.use_signed_output=0;
			t=test(s0,d);
 			u=abs(t.Y-t.X);
			I = find(u < tau);			
			[K,s0.child]=train(s0.child,d); 
			K=K.X;    %% calc kernel				
			w2=s0.alpha'*K*s0.alpha;
			u2=u(I);
			K=K(I,I);
			K2=K .* repmat(r(u2,s,s1,s2),1,length(I));
			[V,D]=eig(K2);
			rho=diag(D);
			n=length(find(rho));
			err=(1-t.Y.*t.X);
			chi=err.*s(err);
			evidence=-lambda*0.5*w2 - sum(chi) - 0.5*sum(lambda+rho) + 0.5*n*log(lambda)			
		elseif strcmp(a.type,'L2')
			if a.use_balanced_C ~= 0
				s0.balanced_ridge=lambda;
			else
				s0.ridge=lambda;
			end;
			[t s0]=train(s0,d);		
			s0.algorithm.use_signed_output=0;
			t=test(s0,d);
			u=abs(t.Y-t.X);
			I = find(u < tau);			
			[K,s0.child]=train(s0.child,d); 
			K=K.X;    %% calc kernel	
			K=add_ridge(K,s0,d);			
			w2=s0.alpha'*K*s0.alpha;					
			K=K(I,I);			
			[V,D]=eig(K);
			rho=diag(D);
			n=length(find(rho));
			chi=0.5*(1-t.Y.*t.X).^2;			
			evidence=-lambda*0.5*w2 - sum(chi) - 0.5*sum(lambda+rho) + 0.5*n*log(lambda)				
		end;	
		gamma=sum(rho./(lambda+rho));	
		if evidence > evidence_max
			evidence_max=evidence;
			lambdamax=lambda;
			gammamax=gamma;			
		end;			
		lambda=gamma/w2; %% gamma_t+1
		if (abs(lambda-lambda_old) < eps) | (i == max_iter)
			finished=1;		
			break;
		end;		
		lambda_old=lambda;
		i=i+1;
	end;
	if strcmp(a.type,'L1')
		a.pbest=1/lambdamax;
	elseif strcmp(a.type,'L2')
		a.pbest=lambdamax;
	end;
	a.posterior=evidence_max/sqrt(gammamax);
else
	while ~finished					
		s0.C=beta/alpha;
		if isinf(s0.C) | (s0.C > 10000)
			s0.C=10000;
		end;
		[t s0]=train(s0,d);
		[K,s0.child]=train(s0.child,d); 
		K=K.X;    %% calc kernel				
		w2=s0.alpha'*K*s0.alpha;			
		t=test(s0,d);
		u=t.Y-t.X;
		I=find(abs(u)-s0.epsilon < eps);
		u2=u(I);
		ri=r(-u2-s0.epsilon,s,s1,s2) + r(u2-s0.epsilon,s,s1,s2);
		K=K(I,I);
		K2=K .* repmat(ri,1,length(I));
		[V,D]=eig(K2);
		rho=diag(D);
		n=length(find(rho));			
		chi=max(0,u-s0.epsilon) + max(0,-u-s0.epsilon);
		ED=sum(chi);
		evidence=-alpha*0.5*w2 - beta*ED - 0.5*sum(log(alpha+beta*rho)) + 0.5*n*log(alpha) + numEx*(log(beta)-log(2)-log(1+s0.epsilon*beta))		
		gamma=sum((beta*rho)./(alpha+beta*rho));
		if evidence > evidence_max
			evidence_max=evidence;
			alphamax=alpha;
			betamax=beta;			
			gammamax=gamma;
			EDmax=ED;
		end;			
		alpha=gamma/w2;		
		b=(2*ED+gamma*s0.epsilon)/(4*s0.epsilon*ED);
		beta=-b + sqrt(b^2 + (2*numEx-gamma)/(2*s0.epsilon*ED));
		if ((abs(alpha-alpha_old) < eps) & (abs(beta-beta_old) < eps)) | (i == max_iter)
			finished=1;		
			break;
		end;		
		alpha_old=alpha;
		beta_old=beta;
		i=i+1;
	end;
	a.pbest=betamax/alphamax;
	dalpha2=2/gamma;
	dbeta2=(1+s0.epsilon*betamax)^2/(betamax*(EDmax*(1+s0.epsilon*betamax)^2+numEx*s0.epsilon));
	a.posterior=evidence_max*sqrt(dalpha2*dbeta2);
end;
