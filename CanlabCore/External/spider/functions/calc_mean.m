
function [r,s]= calc_mean(d,loss_type,as_vector,take_average)

%  CALC_MEAN calculate mean loss across groups of data objects
%   
%  [X,S]=CALC_MEAN(D,L,K) calculates the mean for groups of data objects in D.
%  A loss function L can also be supplied if a loss function has not 
%  already been applied.
%  The (optional) parameter K if set to 0 (default) specifies if
%  the results should be stored in a data object or if set to 1, a vector. 
%  S always returns a vector of standard errors, one for each group.

%  The data objects D should be stored as a cell array of cell arrays
%  such that length(D) is the number of methods to compare and 
%  length(D{i}) is the number of trials, which is equal for all i.
%  If data is not grouped (it is stored as a flat cell array, an attempt
%  is made to group the data automatically (e.g splitting into cv folds).

if nargin<4 take_average=0; end;
if nargin<3 as_vector=0; end;
if nargin<2 loss_type=[]; end;

if isa(d.child{1},'group') 
  n=get_name(d.child{1}.child{1}); 
else  
  n=get_name(d.child{1}); 
end
nocv=0; if isempty(findstr('cv',n)) nocv=1; end;

fast=1; 
if strcmp(loss_type,'confusion_matrix') fast=0; end;


  
%% ------- fast call but doesn't cope with complex objects ------
if ~nocv & fast
 [r,s]=calc_mean_fast(d,loss_type,as_vector,take_average);
 if ~isempty(r) return; end;
end

calc_loss=1; if isempty(loss_type) calc_loss=-1; loss_type='class_loss'; end;
if nargin<2 calc_loss=-1; loss_type='class_loss'; end;

%% -------------------------------------------------------
 
l=length(d.child);             % number of folds
dd=[]; L=[];
for i=1:l
  dd{i}=group2cell(d.child{i}); %% calc elements in each fold
end
m=length(dd{1}); miscell=iscell(dd{1});

%% ---------  take_average_over={0-'guess',1-'inside_leaf',2-'outside_leaf'}


if take_average==0 & length(dd{1})>1 
  if take_average==0   % guess which one
    take_average=1;
    if ~isempty(findstr(get_name(dd{1}{1}),'fold=1')) 
      if ~isempty(findstr(get_name(dd{2}{1}),'fold=1')) 
	take_average=2;
      end
    end
  end
  if take_average==2
    disp('[assuming take_average on outside_leaf]');
  end
end

if take_average==2 & length(dd{1})>1
  dd2=[];
  for i=1:l
    for j=1:length(dd{i})
      dd2{j}{i}=dd{i}{j};
    end
  end
  dd=dd2; dd2=[]; m=l; l=length(dd); 
end


%% -------------------------------------------------------

if calc_loss==-1  %% check to see if we should calc loss
  if miscell name=get_name(dd{1}{1}); else name=get_name(dd{1}); end;
  if isempty(findstr('_loss',name))
    disp('[assuming class_loss]'); calc_loss=1; 
  end
end

%% -------------------------------------------------------

r=[]; 
for i=[1:m]         %% calculate means
  v2=[];
  for k=1:l
    if ~miscell d2=dd{k}; else d2=dd{k}{i}; end;
    if calc_loss==1 d2=loss(d2,loss_type);  end;
    Y = get_y(d2);
    if k==1 v=Y; r{i}=d2;  else v=v+Y; end;
	v2=[v2 Y(1,1)];
  end 
  y=v/l;       % get mean
  s(i)=std(v2)/sqrt(length(v2)); if size(y,2)==1 y(1,2)=s(i); end; % std err
  r{i} = set_y(r{i}, y);
  % ------ remove "fold" word and add new loss value --------- 
  
  y=get_name(r{i});   
  if ~nocv
    t=findstr('fold=',y); y1=y(1:t-1); y2=y(t+4:length(y));
    t=min(findstr('->',y2)); 
    y=[y1 num2str(l) ' folds ' y2(t:length(y2))];
  end  
  t=max(findstr('_loss',y));
  if ~isempty(t)        %% remove loss  function calculation
    t2=max(findstr('>',y(1:t)));
    y=y(1:t2); 
  end
  
  Y = get_y(r{i});
  y=[y  ' ' loss_type '=' num2str(Y(1,1),4) ... 
    '+-' num2str(Y(1,2),4) ...
    ];
  r{i} = set_name(r{i}, y);
end
r=group(r);

%% -------------------------------------------------------

if as_vector
  e=[];
  for i=1:length(r.child)
   Y = get_y(r{i});
   if size(Y)==[1 2]
     e=[e Y(1,1)]; 
   else
     e=[e Y];
   end
  end
  r=e;
end
