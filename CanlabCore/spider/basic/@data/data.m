function data = data(name,X,Y) 

%========================================================================================================  
% Data object
%========================================================================================================
% Stores data into two components X (input) and Y (output).
%
% An object is created with:
%   data('name',X,Y); 
%       or 
%   data(X,Y); (no name)
%       or 
%   data(X); (in case there's no output, e.g clustering)
%
%
% Public attributes:
%  X          -- input matrix        
%  Y          -- output matrix       
%  index      -- original parent indices of examples 
%  findex     -- original parent indices of features 
%
% Public methods:
%  d=get(data,index,featInd)             -- return data,index only (ind)examples,(find)features 
%  x=get_x(data,index,fInd)              -- return x,index only (ind)examples,(find)features 
%  y=get_y(data,index,fInd)              -- return y,index only (ind)examples,(find)features 
%  d=set_x(data,X,indes,featInd)         -- set inputs indexed by (ind)examples,(find)features
%  d=set_y(data,Y,index,featInd)         -- set outputs indxed by (ind)examples,(find)features
%  [indes,featInd]=get_index(d)          -- returns example and feature indices
%  [numEx,vDim,oDim,numCls]=get_dim(d)   -- returns number of examples,features,output dimensions,classes
%
% Import methods:
%  The Data object can import libsvm or arff files for file exchange with weka/libsvm.  
%  Note that arff does not encode what is input or output thus the resulting data object has an empty Y member.
%  Nominal attributes are mapped to an index value and missing attributes are indicated by "NaN".
% 
% Example: 
%  get_x(data([1:5;6:10;11:15]),[1 3],[3:5])
%  dlibsvm=readfrom(data,'libsvm','mylibsvmtrainingfile');
%  darff=readfrom(data,'arff','mywekafile.arff');

%
  if nargin==0
     name='data (empty)';
     data.X=[]; data.Y=[];
     data.index=[]; data.findex=[];
     %data = class(data,'data');
   else 
     data.X=[]; data.Y=[];
     if nargin>=2 data.X=X; end
     if nargin>=3 data.Y=Y; end;
     if ~isa(name,'char') 
       data.Y=data.X;
       data.X=name;
       name=['data'];
     end
     
     X=data.X;
     if isa(X,'cell')
       data.index = [1:length(X)];
       data.findex = [1:length(X{1})];       
     else
       data.index = [1:size(X,1)];
       if isempty(data.index) data.index = [1:size(data.Y,1)]; end;
       data.findex = [1:size(X,2)];       
     end

  end

  p=algorithm(name); p.is_data=1;
  data= class(data,'data',p);
    
  
