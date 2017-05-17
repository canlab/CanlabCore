function a = subsasgn(a,s,value)

  
Q = 'a'; Q2='a';

%% ---- preprocess adding children calls each time there isn't one!!! --
%%  to minimic f.child{1}.child{1}  by just doing f{1}{1} !!

t=s; j=1; isfile=isobject(a); lasttext=0;
for i=1:length(s)

 if ischar(s(i).subs)  %% check if last thing accessed is an object
    isfile = logical(myexist(s(i).subs,'file')==2);
 end

   
  if equal(s(i).type,'{}') | equal(s(i).type,'()')  %% time for child! 
    before=i-1; if i==1 before=1; end; %% look at s(before)
    if ~equal(s(before).subs,'child') & isfile
      % if it exists as file is hopefully an object, isnt child add one!      
      t(j)=s(1); t(j).subs='child'; t(j).type='.';
      j=j+1;
    end
  end
  
  if ischar(s(i).subs) 
    lasttext=j;  
  end
  
  t(j)=s(i); j=j+1;
end
s=t;
%%-------------------------------------------------------------


for i=1:length(s)
  if equal(s(i).type, '.')
    Q=[Q, '.', s(i).subs];
  elseif equal(s(i).type, '()')
    tmp = length(s(i).subs);
    Q=[Q, '('];
    for j=1:(tmp-1)
      Q=[Q, 's(', num2str(i), ').subs{', num2str(j), '}, '];
    end;
    Q=[Q, 's(', num2str(i), ').subs{', num2str(tmp), '})'];
  
  elseif equal(s(i).type, '{}')
    tmp = length(s(i).subs);
    Q=[Q, '{'];
    for j=1:(tmp-1)
      Q=[Q, 's(', num2str(i), ').subs{', num2str(j), '}, '];
    end;
    Q=[Q, 's(', num2str(i), ').subs{', num2str(tmp), '}}'];
  end ;

  if i<lasttext Q2=Q; end;
end;

%% ---- if member doesn't exist search in the children,algorithm.. --------
%%            --- assign value to all possible hits -------


%if ~isfile b=eval(Q2);if isfield(struct(b),name) isfile=1; end; end;
global allow_recursive_subsasgn;
	
if ~isfile & allow_recursive_subsasgn
  name=s(lasttext).subs;  % this is the member we are searching for
  Qend=Q(length(Q2)+length(name)+2:length(Q));

  stack=[]; stack{1}=Q2;
  while length(stack)>0 
        
      cur=stack{1}; stop=0; b=eval(cur);   % eval top line
      if isfield(struct(b),name) 
	Q=[stack{1} '.' name Qend];
  	eval([Q ' = value;']) ;  %% evaluate final answer
	if strcmp(name,'child') stop=1; end; %% shouldnt
                                             %recursively add
                                             %children - it loops
      end; % check for  member

      if ~isfield(struct(b),'alias') & ~stop % as long as not an @algorithm 
	ali=b.algorithm.alias;
	for i=1:2:length(ali)
	  if strcmp(name,ali{i})
	    if isfield(struct(b),ali{i+1}) %% found it!
	      newname=ali{i+1}; 
	      Q=[stack{1} '.' newname Qend];
	      eval([Q ' = value;']) ;  %% evaluate final answer
	    end
	  end
	end
      end
						  
      if isfield(struct(b),'algorithm') & ~stop % check algorithm path 
	stack{length(stack)+1}=[ cur '.algorithm'];
      end
      if isfield(struct(b),'child') & ~stop % check child path
	if ~isempty(b.child)
	  if ~iscell(b.child) 
	    stack{length(stack)+1}=[ cur '.child'];
	  else
	    for i=1:length(b.child)
	      stack{length(stack)+1}=[ cur '.child{' num2str(i) '}'];
	    end
	  end
	end
      end
      if length(stack)>0
	 stack=stack(2:length(stack));
      end
  end
else
  %% simple evaluation.. don't search for multiple versions..
  Q = [Q, ' = value;'];
  eval(Q) ;  %% evaluate final answer
end
			     				     
					     
%%-------------------------------------------------------------





