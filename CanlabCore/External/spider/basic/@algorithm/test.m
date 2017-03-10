function [dat] =  test(algo,dat,lossType)
% 
%               [res]=test(algo,data,loss) 
%   An algorihm algo is trained on the training set data. If a loss
%   function is supplied, it will be calculated.
%   It returns a data object containing the results.
%   If no loss function is supplied, the data object will contain the loss
%   on the training set. If not, it usually contain the label estimates in
%   the X part and the true labels in the Y part of the data object.
%
%   It is also possible to call test(algo), where the emtpy data set is passed 
%   into the algorithm to test it. This is useful for methods that generate their
%   own data.
%   
%   Note: 
%               test(algo) <=> test(algo,[]).

%% Programming note:
% It is used to call testing.m function which doesn't include loss
% function calculations in a child object 

    
if nargin==1 
    dat=[];% <--- data is optional
end; 

if nargin<3 
    lossType=[]; 
end; 

if iscell(dat) 
    dat=group(dat); 
end;
%%<<----test data---->>
if ~isempty(dat) 
    if ~am_i_data(dat)  
        dat=test(dat); %% <--- i.e generate data 
    end
end

%%%% SGE support (note sightly different behaviour from algorithm/train)
if isdeferred(algo)
	if submitted(algo.deferred)
		[algo jobfailed] = waitcollect(algo);
		if jobfailed, return, end
	else
		dat.deferred = qsub(algo.deferred, max(nargout,1), mfilename, algo, dat, lossType);
		return
	end
end
if ~isempty(dat)
if isdeferred(dat)
	[dat jobfailed] = waitcollect(dat);
	if jobfailed, return, end
end 
end
%%%%

e=struct(dat);
%%<<------ if there are multiple datasets as input ---->>
if isa(dat,'group')  & strcmp(e.group,'separate') 
  dat=e.child; % dat seen as @algorithm object => so it is not possible to access child directly
  dat=make_cell(dat);
  res=[];
  for i=1:length(dat)
      [r]=test(algo,dat{i},lossType); 
      res{i}=r;
  end
  dat=group(res);  %% <--- return data/res as set of group
else
    [dat]=testing(algo,dat);
    if ~isempty(lossType) 
      if isa(lossType,'loss')
	dat=train(lossType,dat);
      else
	dat=loss(dat,lossType);
      end
    end;    
    if iscell(dat) 
        dat=group(dat); 
    end;  %%  <--- return data/res as set of group
end



