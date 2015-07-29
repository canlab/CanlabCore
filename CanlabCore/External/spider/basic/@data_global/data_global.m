function data = data_global(name,xx,yy,index,findex) 

%=======================================================================================================  
% DATA_GLOBAL object          
%=======================================================================================================
% Stores data into two components X (input) and Y (output).
% The difference to the data object consists in storing data in global
% variables X and Y to avoid them to be copied on new functions calls. This
% reduces memory overheads.
%
% An object is created with:
%   data_global('name',X,Y);  
%       or 
%   data_global(X,Y); 
%       or
%   data_global(X) (in case there's no name)
%
%   If you want to create it with globals, define X and Y and call
%   data_global('name') or just data_global.
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
% Example: 
%          global X; X=[1:5;6:10;11:15]; get_x(data_global('gi'),[1 3],[3:5])
%
% Note:
% Fast storage of data is achieved by using global X, i.e  X(index,findex).
% xx and yy (optionally) can be stored if different from original. If 
% xx is empty then it is assumed to be equal to global X with given indexes. The same is considered
% for yy ALSO, even if xx instantiated, always global Y(index,:) is used. 
%=========================================================================================================
% Reference : 
% Author    : 
% Link      : 
%=========================================================================================================

  nargins=nargin;
  if nargins==1 & ~ischar(name) xx=name;name='data (empty)';yy=[]; nargins=3; end;
  if nargins==2 yy=[]; nargins=3; end;
  if nargins>=1
    if ~ischar(name)
      if nargins>=2 yy=xx; end; 
      xx=name; name='data (empty)';   
      nargins=3; %nargins+1;
    end
  end
  
  if nargins==3 & (~isempty(xx)  | ~isempty(yy))
       warning(['Creating local copies of data, for globals' ...
		' use data_global without args.']);
     end;

  global X; global Y;
  
  
  if nargins==0
     name='data'; nargins=1;
  end
 
   if nargins==1
     data.myX=[];  %%<--- empty data = X
     data.myY=[];  %%<--- empty data = Y 
   else
     if nargins<2 xx=[]; end;
     if nargins<3 yy=[]; end;
     data.myX=xx; %% <--- data ~= X
     data.myY=yy; %% <--- data ~= Y
   end
   
   if nargins<4       % <--- if there's no indexing
     if isempty(data.myX)
       data.index=[1:size(X,1)];
       if isempty(data.index) data.index = [1:size(data.myY,1)]; end;
       if isempty(data.index) data.index = [1:size(Y,1)]; end;
       data.findex=[1:size(X,2)];
     else
       data.index=[1:size(data.myX,1)];
       if isempty(data.index) data.index = [1:size(data.myY,1)]; end;
       if isempty(data.index) data.index = [1:size(Y,1)]; end;
     end
   else
     data.index=index;
   end
   
   if nargins<5   %% <--- indexes exist, fIndexes do not
     %% do feature indexing
     if isempty(data.myX)
       data.findex=[1:size(X,2)];
     else
       data.findex=[1:size(data.myX,2)];
     end
   else
     data.findex=findex;
   end
      
   
   
  p=algorithm(name); p.is_data=1;
    
  data= class(data,'data_global',p);
   
                                           
   
