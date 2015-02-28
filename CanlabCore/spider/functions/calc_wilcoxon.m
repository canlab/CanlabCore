
function [d]= calc_wilcoxon(d,loss_type,as_vector)

%  CALC_WILCOXON calculate wilcoxon test across data objects
%   
%  [X]=CALC_WILCOXON(D,L,K) 
%  A loss function L can also be supplied if a loss function has not 
%  already been applied.
%  The (optional) parameter K if set to 0 (default) specifies if
%  the results should be stored in a data object or if set to 1, a vector. 

%  The data objects D should be stored as a cell array of cell arrays
%  such that length(D) is the number of methods to compare and 
%  length(D{i}) is the number of trials, which is equal for all i.
%  If data is not grouped (it is stored as a flat cell array, an attempt
%  is made to group the data automatically (e.g splitting into cv folds).



if nargin<3 as_vector=0; end;
if nargin<2 calc_loss=-1; loss_type='class_loss'; else calc_loss=1; end;
if isempty(loss_type) calc_loss=0; end;


if calc_loss==-1  %% check to see if we should calc loss
  dd=group2cell(d.child{1}); 
  if ~(isa(dd,'data') | isa(dd,'data_global') )  dd=dd{1}; end;
  if isempty(findstr('_loss',dd.name))
    disp('[assuming class_loss]'); calc_loss=1; 
  end
end

if calc_loss==1
  d=train(loss(loss_type),d);
end

l=length(d.child);             % number of folds
dd=[];
for i=1:l
  dd{i}=group2cell(d.child{i}); %% calc elements in each fold
end
m=length(dd{1}); miscell=iscell(dd{1});

 
r=[]; xs=[];
for i=1:m
  r1=[];
  for j=1:m
    v=[]; v1=[]; v2=[]; 
    for k=[1:l]
      v1=[v1 dd{k}{i}.Y];
      v2=[v2 dd{k}{j}.Y];
      v=[v dd{k}{i}.Y-dd{k}{j}.Y];
    end 
    %[p side]= wilcoxon_test(v); if side==-1 p=1-p; end; 
    p=signrank(v1,v2); if mean(v1-v2)>0  p=1-p; end;
    r1=[r1 p];
  end
  xs=[xs;r1];
  
   % ------ remove "fold" word and add new loss value --------- 
   
  r{i}=dd{1}{i}.name; y=r{i};
  t=findstr('fold=',y); 
  if ~isempty(y) y1=[y(1:t-1) y(t+4:length(y))]; end;  
  t=max(findstr('->',y));
  if ~isempty(t) %% remove loss  function calculation
    y=[y(1:t-1) '->']; 
  end
  y=[y  ' wilcoxon =['];
  for h=1:m
    y=[y num2str(r1(h)) ' '];
  end; y=[y ']'];
  r{i}=y;
end


d=data(char(r),xs,[]);

if as_vector
  d=d.X;
end

