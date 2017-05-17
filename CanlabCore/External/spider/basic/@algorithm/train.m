function [dat,algor] =  train(algo,dat,lossType)

%
%                  [res,algo] = train(algo, data, loss)
%
%   trains an algorithm using data as training set and (optionally)
%   calculates the loss function. After returning, res contains the results
%   while algo is the updated algorithm, including the learnt model.
%
%   If a loss function is supplied, res contains the loss on the training set.
%   If not, res contains the label estimates in the X part and the true labels 
%   in the Y part of the data object.
%
%   It is also possible to call train(algo), where the empty dataset is
%   passed into the algorithm for training. This is useful for methods,
%   generating their own data (e.g train(chain({toy_data knn})).
% 
%   Note: 
%       train(algo) <=> train(a,[]).
  
%% Programming note:
%  This is used to call training.m function which doesn't include loss
%% function calculations in a child object

    
% INPUT:   single data object (could also be a group of data)
%          single algorithm object (could also be a group)
%
% OUPUTS:  single data object (or a group object of results)
%          single algorithm object (or a group)
   
algo.training_time=cputime;
  
if nargin==1 
    dat=[]; 
end; 

if nargin<3 
    lossType=[]; 
end; 

if iscell(dat) 
    dat=group(dat); 
end;

%<<-----test data----->>
if ~isempty(dat)
    if ~am_i_data(dat)  
        dat=test(dat); %% <----- generate data
    end
end

%%%% SGE support
if isdeferred(algo)
	% note slight difference between this and algorithm/test:
	% we cannot do algo=waitcollect(algo) if submitted, because some spider methods involve calling train on a
	% train algorithm (for example, param/training with use_prev_train set), so there is no way of telling
	% whether train called on a deferred, submitted algorithm is a request to wait
	% for the result of the deferral and then train locally, or to re-train immediately on the cluster.
	def = algo.deferred;
	algo.deferred = [];
	algor = algo;
	def = reset(def, 'description', [mfilename ' ' get_name(algo)]);
	[dat.deferred algor.deferred] = qsub(def, max(nargout, 1), mfilename, algo, dat, lossType);
	return
end
if ~isempty(dat) & isdeferred(dat)
	if isa(algo, 'loss')
		algor = algo;
		dat.deferred = setloss(dat.deferred, algo);
		return
	end
	[dat jobfailed] = waitcollect(dat);
	if jobfailed, algor = algo; return, end
end
%%%%


e=struct(dat); 
% <<-----train multiple data sets----->>
if isa(dat,'group') & strcmp(e.group,'separate')  
    %<<-----dat interpreted as @algorithm object => not accessible yet----->>
    dat=e.child; 
    dat=make_cell(dat);
    algor=[]; 
    res=[];
    for i=1:length(dat)
        if ~(algo.trained==1&algo.do_not_retrain==1),
            [r,al]=train(algo,dat{i},lossType); 
            algor{i}=al;
            res{i}=r;
        else
            res{i} = test(algo,dat{i});  
        end;
    end
    %<<-----data/res is returned as a set of group----->>
    dat=group(res); 
else
  if ~(algo.trained==1&algo.do_not_retrain==1),
	if algo.verbosity==0
		evalc('[dat,algor]=training(algo,dat);');
	else
		[dat,algor]=training(algo,dat);
	end
    algor.trained=1;
    if ~isempty(lossType)  
      if isa(lossType,'loss')
    	dat=train(lossType,dat);
      else
    	dat=loss(dat,lossType);
      end
    end;    
  else
    dat = test(algo,dat); 
    algor=algo; 
  end;
end

%%<<-----correct output type will be determined----->>
if iscell(algor)            
  if length(algor)==1
    algor=algor{1}; % <---single object will be created
  else
    algor=group(algor);% <---create group object
  end
end

algor.training_time=cputime-algor.training_time; %% <--- time needed for training has been measured

