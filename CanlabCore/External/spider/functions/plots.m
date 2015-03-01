
function plots(d,xp,yp,hyper,methods)

% PLOTS  plotting for spider
%   
% plots(d,xp,yp,[hyper],[methods])
%
% Inputs:
%
%  d          -- results in group of data objects
%  xp,yp      -- x-axis and y-axis variables, can be hyperparameters
%                of algorithms, or the type of loss, e.g xp='C', y='class_loss'
%  hyper      -- string of hyperparameters        
%  methods    -- cell array of strings of names methods of interest
%
% Hyperparameters:
%  scatter=0      -- set to 1 for a scatter plot, otherwise tries to draw lines
%  log_scale=0    -- if set to 1 tries to guess when to use a log scale
%  logs=[0 0]     -- set to 1 for x-axis (1st component) or y-axis (2nd)
%                    to be on a log scale
%  no_error_bar=0 -- set to 1 to not draw error bars when given std. dev.
%  ax=[]          -- if nonempty calls axis(ax) to set axis dimensions
%
%  colors=['bkrgymcrgbdkym']          -- default colors to cycle through
%  dots= {'-*','--o','-.',':'};       -- default line styles 
%  scats={'*','o','.','+','x','d','s'};- default scatter plot styles
%
% Examples:
%  r=train(cv(param(svm(kernel('rbf')),'rbf',2.^[-5:5])),gen(toy('l=60')));
%  plots(get_mean(r),'rbf','class_loss','log_scale')    
%  plots(r,'fold','class_loss')    
%
%  [tr a]=train( chain({ param(toy,'l',[10:10:100])  group({svm knn}) }) ) 
%  r=test(a); plots(r,'l','class_loss')
  
if nargin<3 yp='class_loss'; end;
if nargin<4 hyper=[]; end;
if nargin<5 methods=[]; end; 

d=group(group2cell(d));
if (~isempty(findstr(xp,'loss')) | ~isempty(findstr(yp,'loss')) ) 
  if isempty(findstr(d.child{1}.name,'loss'))
    disp('[assuming class_loss]'); 
    d=loss(d,'class_loss');
  end
end


%% ------------------ basic defaults -------------------------
clf
colors=['bkrgymcrgbdkym']; 
dots={'-*','--o','-.',':'};
scats={'*','o','.','+','x','d','s'};
assume=1;                      % make assumptions on methods and losses
scatter=0;                     % scatter plot 
no_error_bars=0;
log_scale=0; logs=[];
xname=xp; yname=yp;  
leg_pos=0; leg_txt=[];
ax=[]; 
%% ------------------- intelligent defaults -------------------
if strcmp(xp,'fold') | strcmp(yp,'fold') assume=0; scatter=1; end;

% labeling and titles
if strcmp(xname,'feat') xname='features'; end;
if strcmp(yname,'feat') yname='features'; end;
xname(xname=='_')=' '; yname(yname=='_')=' ';
%% ------------------------------------------------------------
evalobject=0;
eval_hyper;             %% evalute hyper params
%% ------------------------------------------------------------


%if nargin<2 calc_loss=1; loss_type='default'; else calc_loss=1; end;
%if isempty(loss_type) calc_loss=0; end;

d=d.child; dx=d; dy=d;

[m o t u a ind]=get_methods(d);
if isempty(methods) methods=m; end;
m=[]; for i=1:length(methods) 
  methods{i}(methods{i}=='_')='-';
  m=[m ' ' methods{i}]; end;
disp(['[assuming methods:' m ']']);
if isempty(leg_txt) 
  leg_txt=methods;  
  leg_err_txt=[];
  for i=1:length(methods)
    leg_err_txt{(i-1)*2+1}=[methods{i} ' error bars']; 
    leg_err_txt{(i-1)*2+2}=[methods{i} ]; 
  end
end;

%if assume & ~isempty(strmatch('cv',u)) & strmatch('fold',t{1})
% dx=get_mean(dx); dy=get_mean(dy);
%end

%if length(dx)~=length(dy)  %% problem: grouping occured in one, not the other
%  if length(dx)<length(dy) dy=dx; else dx=dy; end;
%end


m=length(methods); % number of methods
r=[]; index=[];
for i=[1:m]
  x=[];y=[]; xe=[]; ye=[];
  for k=[1:length(dy)]
   if a(k,ind(i))      %% if method i is used in data k this is TRUE

     x1=get_data_value(dx{k},xp);
     y1=get_data_value(dy{k},yp);

     if ~isempty(x1) & ~isempty(y1)
       %% error bars 
       xe1=0;ye1=0;
       if size(x1,2)>1 xe1=x1(1,2); x1=x1(1,1); end;  
       if size(y1,2)>1 ye1=y1(1,2); y1=y1(1,1); end;  
       xe=[xe xe1]; ye=[ye ye1];
       x=[x x1]; y=[y y1];
     end
   end
  end 
  if i==1  
    hold off; 
    if log_scale & isempty(logs)   %% log scale stuff
       if max(x)>10 logs(1)=1; end;
       if max(y)>10 logs(2)=1; end;
    end
    if length(logs)<2 logs(2)=0; end;
  end;
  if logs(1) x=log(x); end; if logs(2) y=log(y); end;

  c=colors(mod(i-1,length(colors))+1);
  if ~scatter
   d=dots{mod(i-1,length(dots))+1};
  else
   d=scats{mod(i-1,length(scats))+1};
  end
  if sum(ye)==0 | no_error_bars
    plot(x,y,[c d ]); hold on;
  else
    errorbar2(x,y,ye,-ye,[c d ]); hold on;
    leg_txt=leg_err_txt;
  end
end

if logs(1) xname=['log ' xname]; end;
if logs(2) yname=['log ' yname]; end;
xlabel(xname); ylabel(yname);

legend(char(leg_txt),leg_pos);
if ~isempty(ax) axis(ax); end;

 
%% ---------------------------------------------------------------------


function val=get_data_value(d,v)
  
vv=length(v);
if  vv>4 & ~isempty(strcmp('_loss',v(vv-5:vv)))
  val=d.Y;
else
  v=[v '=']; n=d.name;
  f=max(findstr(n,v))+length(v);
  if isempty(f) val=[]; return; end;
  p=[]; 
  while 1
    if f>length(n) break; end;
    if n(f)~='=' & n(f)~=' '
      p=[p n(f)];
    else
      break;
    end
    f=f+1;
  end
  val=str2num(p);
end

