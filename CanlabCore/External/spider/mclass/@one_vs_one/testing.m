function dat =  testing(a,dat)

%disp(['testing ' get_name(a) '.... '])

[numEx,vDim] = get_dim(dat);
oDim=a.nrofclasses;

 datTemp=dat; 
 childIndex=0;
  for i=1:oDim,
    for j=(i+1):oDim,
        childIndex=childIndex+1;
        datTemp=set_name(datTemp,['Machine ' num2str(childIndex)]);
        yEst(:,childIndex) = get_x(test(a.child{childIndex},datTemp));    
    end;
  end;
%% compute the score
childIndex=0;
score = zeros(numEx,oDim);
for i=1:oDim,
    for j=(i+1):oDim,
        childIndex=childIndex+1;
        score(:,i) = score(:,i) + sign(yEst(:,childIndex)); 
        score(:,j) = score(:,j) - sign(yEst(:,childIndex));
    end;
end;
%% Take the max:
[dummy, ySort] = sort(-score,2);
yMax = ySort(:,1);
%% If there is some ties, we smooth the output of each classifier
%% by a hyperbolic tangeant and compute the real score for each classes
temp = find(dummy(:,1)==dummy(:,2));
if ~isempty(temp),
  childIndex=0;
  score = zeros(length(temp),oDim);
  for i=1:oDim,
    for j=(i+1):oDim,
        childIndex=childIndex+1;
        score(:,i) = score(:,i) + tanh(yEst(temp,childIndex)); 
        score(:,j) = score(:,j) - tanh(yEst(temp,childIndex));
    end;
  end;
  [dummy, ySort] = sort(-score,2);
  yMaxTemp = ySort(:,1);
  yMax(temp) = yMaxTemp;
end;
yEst = -ones(numEx,oDim);
 for i=1:length(yMax),
    yEst(i,yMax(i))=1;
 end;
 
dat=set_x(dat,yEst); 
dat=set_name(dat,[get_name(dat) ' -> ' get_name(a)]); 

