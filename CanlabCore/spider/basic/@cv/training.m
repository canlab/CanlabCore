function [retRes,retAlgo] = training(retAlgo,dat,lossType)
  
  
%  Does cross validation and returns the results as a data object
%  plus the updated group (after training models).
  
  rndState = rand( 'state');
  if retAlgo.algorithm.verbosity
    disp(['training ' get_name(retAlgo) '.... '])
  end
  if nargin==1 
    dat=data; % <-- in case data is not specified 
  end;          
  
  if nargin<3 
    lossType=[];  % <-- in case loss is not specified
  end 
  
  untrained = retAlgo.child{1};
  if isa(untrained,'cell') 
    untrained = group(untrained); % <-- convert retAlgo.child{1} to group
  end;
  
  [num,vecDim,outDim]=get_dim(dat);
  orig_name=get_name(dat);
  retRes = {};
  
  retAlgo.trialbytrial=[];
  for loop=1:retAlgo.repeats % <-- do cross validation vecDim times

    rand('seed',loop);                
    if retAlgo.balanced     % <-- balance the number of positives/negatives labels in folds
      labels=get_y(dat); 
      fin1=find(labels==1); 
      fin2=find(labels==-1);
      
      pe1=randperm(sum(labels==1)); 
      pe1=fin1(pe1);
      
      pe2=randperm(sum(labels==-1)); 
      pe2=fin2(pe2);
      perm=[pe1 ;pe2];   
    else
      perm=randperm(num);
    end 

    trialbytrial = [];
    for i=1:retAlgo.folds
      if retAlgo.train_on_fold==0
        tst=perm(i:retAlgo.folds:num);
        trn=setdiff(1:num,tst);
      else  % <--- take fraction as training data
        trn=perm(i:retAlgo.folds:num);
        tst=setdiff(1:num,trn);
      end
      
      if retAlgo.store_all 
        indx=i+(loop-1)*(retAlgo.folds); 
      else 
        indx=1; 
      end;
      
      str=[' fold=' num2str(i+(loop-1)*(retAlgo.folds))]; 
      if retAlgo.output_train_error==0  % <--- testing on left out fold
        dat=set_name(dat,[ orig_name ' -> cv' str]);      
      else
        dat=set_name(dat,[ orig_name ' -> cv' str ' output_train_error=1']);   
      end

      if isfield(struct(untrained), 'test_on_the_fly')
        if isstruct(untrained.test_on_the_fly)
          %% untrained.test_on_the_fly.data = set_name(get(dat, tst), [orig_name ' -> cv' str]);
          % line above fails because subsasgn has suddenly stopped working recursively
          % work around it as follows:
          sss = untrained.test_on_the_fly;
          sss.data = set_name(get(dat, tst), [orig_name ' -> cv' str]);
          untrained.test_on_the_fly = sss;
          % some algorithms train another algorithm multiple times in the course of
          % their own training: for example, RFE trains a classifier once per
          % feature elimination step. It is often handy to evaluate cv-error at each
          % stage: previously this involved either re-training the individual classifiers
          % afterwards (time-consuming), or storing a whole series of trained
          % classifiers (memory-consuming). This solution gives the wrapper
          % a sneak preview of the test data for this fold, so that cv error can be
          % evaluated as it goes. Of course, the sneak preview should NOT be used
          % for training!
        end
      end
      if retAlgo.algorithm.verbosity
        str = sprintf('training cv fold %d of %d', i, retAlgo.folds);
        if retAlgo.repeats > 1, str = sprintf('%s (repeat %d of %d)', str, loop, retAlgo.repeats); end
        disp(str)
      end
      
      if retAlgo.output_train_error
        [res, retAlgo.child{indx}] = train(untrained, get(dat, trn), lossType);
      else
        [res, retAlgo.child{indx}] = traintest(untrained, dat, trn, tst, lossType);
      end
      if retAlgo.store_trialbytrial
        trialbytrial = [trialbytrial; tst(:) repmat(i, length(tst), 1) res.Y res.X];
      end
      %	[res retAlgo.child{indx}]=train(untrained,get(dat,trn));  
      %    if retAlgo.output_train_error==0  % <--- test on left out fold
      %      [res]=test(retAlgo.child{indx},get(dat,tst)); 
      %    end            
      %    if ~isempty(lossType)        % <--- calculate loss on results
      %      res=loss(res,lossType,[],1); 
      %    end
      
      if ~isa(res, 'cell'), res = {res}; end
      retRes = [retRes res(:)'];	
      
    end
    if retAlgo.store_trialbytrial
      retAlgo.trialbytrial = cat(3, retAlgo.trialbytrial, sortrows(trialbytrial));
    end
  end

  retRes=group(retRes); 

  rand( 'state', rndState);  % restore the random generator