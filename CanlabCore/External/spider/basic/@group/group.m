
function a = group(in,hyper) 

%===========================================================================   
% Group object
%===========================================================================  
%  A=GROUP(I,[G]) returns a group object initialized with 
%  a cell array of algorithms I and an optional grouping type G.  
% 
%  This is used to collect a group of algorithms.
%
%  Examples:  group({svm knn c45}),
%             get_mean(train(cv(group({svm knn })),train(toy)))
%
% 
%  Hyperparameter
%  group='all'           This parameter stores grouping type:
%                        'one_for_each', 'all' or 'separate' (default='all').
%                        When results of training and testing are
%                        output by a (grouped) algorithm,
%                        they are stored in the same group type given
%                        here. This is important for objects like chain, 
%                        and get_mean which do not deal with data objects
%                        separately.
%               
%
%      GROUP TYPE            
%      'all':               each item is passed independently and separately
%                           to all training objects, e.g
%                             d=gen(toy); a=group({knn svm});
%                             train(a,group({d d},'separate'))
%                           gives 4 outputs.
%
%      'one_for_each':      the n^th data item is passed to the
%                           n^th training object in an group object, e.g
%                             d=gen(toy);g=group({d d},'one_for_each');
%                             a=group({knn svm}); train(a,g)
%                           gives 2 outputs.
%
%       'all':              passes all the data objects to a
%                           single training object, e.g 
%                             d=gen(toy); a=group({knn svm}); r=train(a,d);
%                             train(get_mean,group(r,'all'))
%                           gives 1 ouput.
%       
% [NOTE: It is also possible to use the transpose operator with
% group, e.g: r=group({ {svm svm} {knn knn} })' will give
% the same result as group({ {svm knn} {svm knn}  })]  
%===========================================================================
% Reference : 
% Author    :   
% Link      : 
%===========================================================================
  a.child=[];
  a.group='separate';
  if nargin==2
    a.group=hyper;
    if isa(in,'group') 
      a=in; a.group=hyper; return;
    end
  end

  
  if nargin>0
    if length(in)==1&~iscell(in), a.child{1}=in; else a.child=in; end;
  end
   
  datacount=0;
  for i=1:length(a.child) 
    if isa(a.child{i},'cell') a.child{i}=group(a.child{i}); end;
%    if a.child{i}.algorithm.is_data datacount=datacount+1; end;
    if am_i_data(a.child{i}) datacount=datacount+1; end;
  end

  if datacount==length(a.child) 
    name='data/results';   p=algorithm(name); p.is_data=1; 
  else 
    name='group';   p=algorithm(name); 
  end;  

    a= class(a,'group',p);


