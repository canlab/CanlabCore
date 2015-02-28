function results =  test_slow(algo,dt)

%
%
%            results =  test(alg,data)
%
% Generates labels using knn.
  
  if isa(dt,'data_global') %% use data_global in case of low memory
    results=test_fast(algo,dt);
  else
    
    [x  y]  = get_xy(algo.dat);
    [xTemp yTemp] = get_xy(dt);

    if ~strcmp(algo.ker{1},'linear'),
      kerMaTemp=get_kernel(dt,x,algo.ker)';
      zeroKerMa = zeros(size(kerMaTemp));
      
      for i=1:size(kerMaTemp,1),
        kerTemp = get_kernel(dt,xTemp(i,:),algo.ker);
	    zeroKerMa(i,:) =  kerTemp(i); 
      end; 
      
      zeroKerMa2 = zeros(size(kerMaTemp));
      for i=1:size(kerMaTemp,2),
	    kerTemp = get_kernel(algo.dat,x(i,:),algo.ker);
	    zeroKerMa2(:,i) =  kerTemp(i); 
      end;
      
      kerMa = (zeroKerMa + zeroKerMa2 - 2*kerMaTemp)';
    else
      kerMa=dist(algo,get_x(algo.dat),xTemp');
    end;

    if size(y,2)==1,
      for i=1:size(kerMa,2) 
    	[val,valInd]=sort(kerMa(:,i)); 
        if length(find((y==1)|(y==-1)))~=size(y,2),
            yEst(i) = mean(y(valInd(1:algo.k)));    
        else,
        	tabul=tabulate(y(valInd(1:algo.k))+2); %% <--- calculation of most frequently occuring label 
	        [maxComp maxCompInd]=max(tabul(:,2)); 
	        yEst(i)=tabul(maxCompInd,1)-2;
        end;
      end;
      yEst=yEst';
    else
      yEst = -ones(size(kerMa,2),size(y,2));
      for i=1:size(kerMa,2) 
    	[val,valInd]=sort(kerMa(:,i)); 
    	tabul=sum(y(valInd(1:algo.k),:),1); 
    	[maxComp maxCompInd]=max(tabul); 
    	yEst(i,maxCompInd)=1;
      end;
    end;
        
    
    results=data([get_name(dt) ' -> ' get_name(algo)],yEst,yTemp);
    
    
  end
