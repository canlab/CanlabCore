
%% used to set hyperparameters -- needed for every object! 

if  isempty(hyper) | (iscell(hyper) & isempty(hyper{1}) & length(hyper)==1 ) return; end;

% ------- try to split up input via semi-colons ------
if ~iscell(hyper) hyper={hyper}; end; % make it a cell
extra=[];
for j=1:length(hyper)
 if ischar(hyper{j})
   h=hyper{j}; ex=[];
   if h(length(h))==';' h=h(1:length(h)-1); end;
   f=find(h==';');  f=[0 f length(h)+1 ];
   for i=1:length(f)-1
     ex{i}=h(f(i)+1:f(i+1)-1);
   end
     hyper{j}=ex{1};
     a1=extra; if ~isempty(a1) a1=make_cell(a1); end;
     a2=ex(2:length(ex)); if ~isempty(a2) a2=make_cell(a2); end;
     extra=[a1 a2];
 end
end
if ~isempty(extra)
  hyper=[hyper extra];
end
for j=1:length(hyper)
     if iscell(hyper{j}) hyper{j}=hyper{j}{1}; end;
%% --------------- add '=1' if no equals -------------
     if ischar(hyper{j})                
      if isempty(findstr('=',hyper{j}))
       hyper{j}=[hyper{j} '=1'];
      end
     end
%% --------------- add algorithm calls -------------------------
%% -------------------------------------------------------------
end

% ------------if only one hyper no cell ------------
if length(hyper)==1& (isa(hyper,'kernel') | isa(hyper,'distance'))
      a.child=hyper;
else
% ------------ evaluate input --------------

 for i=1:size(hyper,2)
      if isa(hyper{i},'kernel') | isa(hyper{i},'distance')
	  a.child=hyper{i};
      else
           value=hyper{i};
           if ischar(value)
             value(find(value=='"'))=char(39);
           end
	   if exist('evalobject') & evalobject==0
	     eval([value ';']);
	   else
         eval(['a.' value ';']);
	   end
      end
     end;   
    
end;
