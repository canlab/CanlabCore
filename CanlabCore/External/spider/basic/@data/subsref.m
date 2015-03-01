function b = subsref(a,s)
%    B = SUBSREF(A,S) is called for the syntax A(I), A{I}, or A.I
%    when A is an object.  S is a structure array with the fields:
%        type -- string containing '()', '{}', or '.' specifying the
%                subscript type.
%        subs -- Cell array or string containing the actual subscripts.
%

% File:        @data/subsref.m
%
% Author:      Gunnar R"atsch, Alex Smola, Jason Weston
% Created:     02/11/98
% Updated:     02/21/98 05/08/00 07/12/02
%
% This code is released under the GNU Public License
% Copyright by GMD FIRST and The Australian National University

Q = 'a'; Q2='a';

%% ---- preprocess adding children calls each time there isn't one!!! --
%%  to minimic f.child{1}.child{1}  by just doing f{1}{1} !!

t=s; j=1; isfile=isobject(a); lasttext=0;
for i=1:length(s)

  if ischar(s(i).subs)  %% check if last thing accessed is an object
    isfile=logical(myexist(s(i).subs,'file')==2);
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


%% This is the core of the code to convert matlab stuff to eval function

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




%% ------------------ @ data special bit --------------

if ~isfile   %% check for special data members:  X and Y
  name=s(lasttext).subs;  % this is the member we are searching for  
  if strcmp(name,'X') | strcmp(name,'Y')
    if equal(s(length(s)).type,'()') %% indexes into X
      i=s(length(s)).subs{1}; fi=s(length(s)).subs{2}; 
      Q=['get_x(' Q2 ',i,fi);'];
    else
      Q=['get_x(' Q2 ');'];
    end
    if strcmp(name,'Y') Q(5)='y'; end; %% change to get y instead!
    b=eval(Q);return;  %% evaluate final answer
  end
end

%% ------------------ END @ data special bit --------------





%% ---- if member doesn't exist search in the children,algorithm.. --------
%% -- differs from subsasgn in that it only returns first one   ----------

global recursive_subsasgn_off;
if isempty(recursive_subsasgn_off) recursive_subsasgn_off=0;end;
	
if ~isfile & (~recursive_subsasgn_off)
  name=s(lasttext).subs;  % this is the member we are searching for
  Qend=Q(length(Q2)+length(name)+2:length(Q));

  notfound=1;           
  stack=[]; stack{1}=Q2; 
  while notfound %% keep searching into children..

      cur=stack{1};
      b=eval(cur);   % eval top line
      if isfield(struct(b),name) notfound=0; end; % check for  member
						  
      if notfound                      %% check for aliases 
	if ~isfield(struct(b),'alias') % as long asn't currently in @algorithm 
	  ali=b.algorithm.alias;
	  for i=1:2:length(ali)
	    if strcmp(name,ali{i})
	      if isfield(struct(b),ali{i+1}) %% found it!
		name=ali{i+1}; notfound=0;
	      end
	    end
	  end
	end
      end
      
      if notfound
	if isfield(struct(b),'algorithm') % check algorithm path 
	  stack{length(stack)+1}=[ cur '.algorithm'];
	end
	if isfield(struct(b),'child') % check child path
	  if ~iscell(b.child)
	    stack{length(stack)+1}=[ cur '.child'];
	  else
	    for i=1:length(b.child)
	      stack{length(stack)+1}=[ cur '.child{' num2str(i) '}'];
	    end
	  end
	end
	
	stack=stack(2:length(stack));   % remove top of stack
      end
  end
  Q=[stack{1} '.' name Qend];
  
end

%%-------------------------------------------------------------



b=eval(Q) ;  %% evaluate final answer

