function [results,a] =  training(a,d)

% [results,algorithm] =  training(algorithm,data,loss)

% set various parameters
R=a.dist;
Rd = a.distd;
N=get_dim(d);
tolfun = a.tolfun;
tol = a.tol;

if ~strcmp(a.mode,'custom')
  % caluculate pairwise distance
  D=calc(a.child,d);
  K=(D<=R);
  Kd = []; Y = [];
end

if strcmp(a.mode,'discriminative') || strcmp(a.mode,'shared')
  Kd = (D<=Rd); 
  Y= get_y(d); 
end


%
% Just optimize once. Use L1 approximation to LVQ
%
if strcmp(a.optmethod,'l1lvq')
    % perform standard LVQ
    if strcmp(a.mode,'standard')
        % calculate gamma
        gamma=sum(K)';
        % gamma=ones(N,1);
        A=K;
        b=[ones(N,1)];
        [X,fval]=optimize(a,gamma,-A,-b,zeros(size(A,2),1),inf*ones(size(A,2),1));

    elseif strcmp(a.mode,'custom')
        % calculate gamma
        A = get_x(a.child);
	% b was a additive
	% constant in the OR, so it has to be
	% substracted from 1
	b = ones(size(A,2),1) - get_y(a.child); 	
	gamma = ones(size(A,2),1);
        [X,fval]=optimize(a,gamma,-A,-b);
    elseif strcmp(a.mode,'discriminative')
        % calculate gamma
        gamma=sum(K)';
        % gamma=ones(N,1);
        A=K;
        A(Y==1,:) = A(Y==1,:).*repmat(double(Y==1)',length(find(Y==1)),1);
        A(Y==-1,:) = A(Y==-1,:).*repmat(double(Y==-1)',length(find(Y==-1)),1);
        b = ones(N,1);
        b(Y == 1,:) = b(Y == 1,:) - Kd(Y==1,:)*double(Y==-1);
        b(Y == -1,:) = b(Y == -1,:) - Kd(Y==-1,:)*double(Y==1);

        [X,fval]=optimize(a,gamma,-A,-b,zeros(size(A,2),1),inf*ones(size(A,2),1));
    elseif strcmp(a.mode,'shared')
        % calculate gamma
        gamma=sum(K)';
        % gamma=ones(N,1);
        A=K;
        A(Y==1,:) = A(Y==1,:).*repmat(double(Y==1)',length(find(Y==1)),1);
        A(Y==-1,:) = A(Y==-1,:).*repmat(double(Y==-1)',length(find(Y==-1)),1);
        b = ones(N,1);
        b(Y == 1,:) = b(Y == 1,:) - prod(ones(length(find(Y==-1)),length(find(Y==1)))-Kd(Y==1,Y==-1)')';
        b(Y == -1,:) = b(Y == -1,:) - prod(ones(length(find(Y==1)),length(find(Y==-1)))-Kd(Y==-1,Y==1)')';

        [X,fval]=optimize(a,gamma,-A,-b,zeros(size(A,2),1),inf*ones(size(A,2),1));
    end
end






%
% Optimize several times. Use L0 approximation to LVQ.
%
if strcmp(a.mode,'discriminative') || strcmp(a.mode,'shared')
  K = double(K);
end



if strcmp(a.optmethod,'l0arom')
  % some initializations
  z = ones(1,N); % downweighting factors
  f = ones(N,1); % objective function
  b = ones(N,1); % >= 1 for w minimization
  sp = 1;   spold = 0; % sparsity trackers
  
  if strcmp(a.mode,'discriminative')
    b(Y == 1,:) = b(Y == 1,:) - Kd(Y==1,:)*double(Y==-1);
    b(Y == -1,:) = b(Y == -1,:) - Kd(Y==-1,:)*double(Y==1);
  elseif strcmp(a.mode,'shared')
    b(Y == 1,:) = b(Y == 1,:) - prod(ones(length(find(Y==-1)),length(find(Y==1)))-Kd(Y==1,Y==-1)')';
    b(Y == -1,:) = b(Y == -1,:) - prod(ones(length(find(Y==1)),length(find(Y==-1)))-Kd(Y==-1,Y==1)')';
  elseif strcmp(a.mode,'custom')
    % calculate gamma
    K = get_x(a.child);
    b =  ones(size(K,1),1) - get_y(a.child);%=b+Y 	
  end

  if strcmp(a.mode,'discriminative') || strcmp(a.mode,'shared')
    K(Y==1,:) = K(Y==1,:).*repmat(double(Y==1)',length(find(Y==1)),1);
    K(Y==-1,:) = K(Y==-1,:).*repmat(double(Y==-1)',length(find(Y==-1)),1);
  end


  A = K;
  lb =  zeros(size(A,2),1); % lower bound
  ub= inf*ones(size(A,2),1); % upper bound


  X = z';
  while abs(spold-sp) > 0
    A = K.*repmat(z,size(K,1),1);
    disp(['Sparsity level is ' num2str(sp)]);
    [X,fval]=optimize(a,f,-A,-b,lb,ub);     
    z = z.*X(1:N)';

    % ------- check if sparsity level has changed -------
    X = X(1:N); [maxS,I]=sort(-abs(X));
    maxV=max(X); minV=min(X);
    maxS=maxS/maxV;
    spold = sp;
    sp = length(find(abs(maxS) > a.cutoff))/N;
  end
  X = X(1:N);
end



% -----------------------------------------------------------------------
% Process outputs, store and test
% format outputs
[maxS,I]=sort(-abs(X));

maxV=max(X);
minV=min(X);
maxS=maxS/maxV;


C=find(abs(maxS) > a.cutoff);
a.keep=get(d,I(C));
a.alpha=X(I(C));
if a.test_on_trainingset == 1 && ~strcmp(a.mode,'custom')
    results=test(a,d);
end
results = d;



function is = feasible(A,x,b)
    x(x>0) = 1;
    if A*x <= b
       is = 1;
    else
        is = 0;
    end
