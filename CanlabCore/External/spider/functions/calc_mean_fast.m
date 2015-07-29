
function [r,s]= calc_mean_fast(d,loss_type,as_vector,take_average)

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


if nargin<3 as_vector=0; end;
if nargin<2 loss_type=[]; end;
calc_loss=1;
if isempty(loss_type) calc_loss=-1; loss_type='class_loss'; end;

%% -------------------------------------------------------

miscell=0; if isa(d.child{1},'group') miscell=1; end; 

if calc_loss==-1  %% check to see if we should calc loss
  if miscell name=d{1}{1}.name; else name=d{1}.name; end;
  if isempty(findstr('_loss',name))
    disp('[assuming class_loss]'); calc_loss=1; 
  end
end

L=[]; r=[];
if isa(d.child{1},'group')
  r=[]; l=length(d.child);
  for i=1:l
    d2=d.child{i}.child;
    for j=1:length(d2)
      if i==1 r{j}=d2{j}; end;
      if isa(d2{j},'group') 
	r=[];s=[]; 
	disp(['[cant do fast get_mean: complex objects]']); 
	return; 
      end;
      if calc_loss==1 
	L(i,j)=get_y(loss(d2{j},loss_type));
      else
	L(i,j)=d2{j}.Y;
      end;
    end; 
  end
else
  r{1}=d.child{1}; l=length(d.child);
  for i=1:l
    if calc_loss==1 
      L(i,1)=get_y(loss(d.child{i},loss_type));
    else
      L(i,1)=d.child{i}.Y;
    end;
  end
end


%% ---------  take_average_over={0-'guess',1-'inside_leaf',2-'outside_leaf'}

if nargin<4 take_average=0; end;
if take_average==0 & miscell
  if take_average==0   % guess which one
    take_average=1;
    if ~isempty(findstr(d{1}{1}.name,'fold=1')) 
      if ~isempty(findstr(d{2}{1}.name,'fold=1')) 
	take_average=2;
      end
    end
  end
  if take_average==2
    disp('[assuming take_average on outside_leaf]');
  end
end

if take_average==2 & miscell
  L=L';
  r=[];
  for i=1:l
    r{i}=d.child{i}.child{1};
  end
end



%% -------------------------------------------------------

m=mean(L); s=std(L)/sqrt(size(L,1));
for i=1:length(r)
  y=r{i}.algorithm.name;
  t=findstr('fold=',y); y1=y(1:t-1); y2=y(t+4:length(y));
  t=min(findstr('->',y2));  y=[y1 num2str(size(L,1)) ' folds ' y2(t:length(y2))];
  t=max(findstr('_loss',y));
  if ~isempty(t)        %% remove loss  function calculation
    t2=max(findstr('>',y(1:t)));
    y=y(1:t2-3); 
  end
 r{i}.algorithm.name=[y ' -> ' loss_type '=' num2str(m(i)) '+-' num2str(s(i))];
 r{i}.Y=[m(i) s(i)]; r{i}.X=[];
end 
r=group(r);


if as_vector
  e=[];
  for i=1:length(r.child)
   if size(r{i}.Y)==[1 2]
     e=[e r{i}.Y(1,1)]; 
   else
     e=[e r{i}.Y];
   end
  end
  r=e;
end


