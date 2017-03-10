
function [r]= calc_group(d,param)
  
% r=GROUP(data,param)
%
% Group data objects into cell arrays : the idea is that one groups
% data objects into a cell array of cell arrays, so that each cell
% array is a group with several members. Then various functions can be
% applied to this grouping e.g the wilcoxon test or calculating the mean
% error of each group.
%
% Data can be grouped in the follow ways (according to the value of param):
% 'cv'     : groups according to algorithms using cv folds 
% 'rotate' : rotates cell array so groups become members and members are groups
% 'flat'   : flatten cell array of cell arrays into single array
% only groups objects according to algorithms used in cv

%% defaults
if nargin==1
 if length(d{1})==1  %% flat structure -- assume cv
  k=findstr(d{1}.name,'fold='); 
  if ~isempty(k)   
    param='cv';
    disp('[group: assuming cv folds]');
  else
    param='flat';
    disp('[group: assume single set]');
  end
 else                %% non-flat -- assume rotation of tree
  param='rotate';
  disp('[group: assuming rotation of members]');
 end
end



if isa(d,'group') %% first extract all algorithms from inside it
  d=group2cell(d);
end

if length(d)==1
 r=d; return; % one object -- nothing to do
end

if strcmp(param,'cv')
count=0;
for i=1:length(d)
 
     n=d{i}.name;        %%------find correct place to put data
     k=findstr(n,'fold=');
     if isempty(k) 
       place=1;
     else
       j=5; c=[];
       while 1
         if isempty(str2num(n(k+j))) break; end;
         c=[c n(k+j)]; j=j+1;
       end
       place=str2num(c);
     end
    if length(count)<place
       count(place)=1;
     else
       count(place)=count(place)+1;
    end
    place=count(place);
    pos(i)=place;
end
end
if strcmp(param,'rotate')
d2=[]; k=1;
   for i=1:length(d)
     for j=1:length(d{i})
       d2{k}=d{i}{j}; 
       pos(k)=j; k=k+1;
     end
   end 
   d=d2;
end
if strcmp(param,'flatten') | strcmp(param,'flat')
  if length(d{1})==1 r=d; return; end;
  d2=[]; k=1;
   for i=1:length(d)
     for j=1:length(d{i})
       d2{k}=d{i}{j}; 
       pos(k)=k; k=k+1;
     end
   end 
   d=d2;
end
r=[];
for i=1:length(d)        %% assign data into correct place
     while 1             
      if isempty(r)
        r{pos(i)}=d{i}; break; % make_cell(d{i}); break; 
      end
      if length(r)<pos(i)
        r{pos(i)}=d{i}; break; %make_cell(d{i}); break; 
      end
      if length(r{pos(i)})==1 r{pos(i)}=make_cell(r{pos(i)}); end; 
      r{pos(i)}={r{pos(i)}{1:length(r{pos(i)})} d{i} }; break;
     end
end
